from pathlib import Path
import json,csv,time,argparse,itertools,concurrent.futures
import numpy as np
from final_model import ROOT,CONTRACTS,Planner,design,bm
def evaluate(job):
    base,output=job
    out=Path(output)/('fleet_'+base['id']+'.json')
    if out.exists():return base['id'],'cached'
    spec=design();p=bm.Physics(base);rows=[];costs=[];observations={};checks=0
    for channel in spec['channels']:
        m=Planner(p,channel);m.prepare_noinfo()
        observations[channel]={c:dict(cells=o.cells,likelihood=[x.tolist() for x in o.likelihood],masks=[x.tolist() for x in o.masks]) for c,o in m.channels.items()};checks+=sum(o.corner_checks for o in m.channels.values())
        for c in CONTRACTS:
            m.configure([c],equal=spec['equal_cost'])
            for method in ['exact','myopic']:
                rows.append(dict(base=base['id'],group=base['group'],channel=channel,contract=c,method=method,te=[ss[0]['te'] for ss in base['states']],**m.run(method,trace=method=='exact')))
        if base['rep']==0:
            for sample,dur in itertools.product(spec['cost_grid']['sample'],spec['cost_grid']['duration']):
                prices=(spec['cost_grid']['setup'],sample,dur)
                for c in CONTRACTS[:-1]:
                    m.configure([c],prices=prices);costs.append(dict(base=base['id'],group=base['group'],channel=channel,sample=sample,duration=dur,contract=c,method='exact',**m.run('exact')))
                m.configure(CONTRACTS[:-1],prices=prices)
                for method in ['exact','myopic']:costs.append(dict(base=base['id'],group=base['group'],channel=channel,sample=sample,duration=dur,contract='joint',method=method,**m.run(method)))
    out.write_text(json.dumps(dict(base=base['id'],rows=rows,costs=costs,observations=observations,corner_checks=checks,physics_seconds=p.seconds),indent=2));return base['id'],len(rows)+len(costs)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2);ap.add_argument('--limit',type=int);ap.add_argument('--output',type=Path,default=ROOT/'results');args=ap.parse_args();bases=design()['bases'][:args.limit] if args.limit else design()['bases']
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(evaluate,[(b,str(output)) for b in bases]):print(result,flush=True)
    for typ in ['rows','costs']:
        rows=[{k:v for k,v in r.items() if k!='leaves'} for b in bases for r in json.loads((output/('fleet_'+b['id']+'.json')).read_text())[typ]]
        with (output/('targeted_'+typ+'.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=sorted(set(k for r in rows for k in r)));w.writeheader();w.writerows(rows)
    print('DONE',len(bases),'fleets',flush=True)
if __name__=='__main__':main()
