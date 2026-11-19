import json,csv,collections,itertools,math
import numpy as np
from final_model import ROOT,design,Planner,bm,CONTRACTS
def stat(x):
    a=np.array(x,dtype=float);return dict(n=len(a),mean=float(a.mean()),min=float(a.min()),max=float(a.max()),positive=int(sum(a>1e-7)),zero=int(sum(abs(a)<=1e-7)),negative=int(sum(a< -1e-7)))


def closure_difference(first, second, equal_support=False):
    if equal_support:
        return 0.0, 0.0
    lower = math.nextafter(first['closure_lower'] - second['closure_upper'], -math.inf)
    upper = math.nextafter(first['closure_upper'] - second['closure_lower'], math.inf)
    # The second support is a subset of the first for both comparisons below.
    assert upper >= 0
    return max(0.0, lower), upper


def main():
    docs=[json.loads(p.read_text()) for p in sorted((ROOT/'results').glob('fleet_*.json'))];rows=[r for d in docs for r in d['rows']];costs=[r for d in docs for r in d['costs']];summary=dict(groups={},counts=dict(fleets=len(docs),matched_policy_rows=len(rows),cost_policy_rows=len(costs),all_boundary_corner_checks=sum(x['corner_checks'] for x in docs)))
    groups=sorted(set(r['group'] for r in rows))
    for group in groups:
        summary['groups'][group]={}
        for ch in ['ideal','bounded']:
            rr=[r for r in rows if r['group']==group and r['channel']==ch and r['method']=='exact'];by={(r['base'],r['contract']):r for r in rr};bs=sorted(set(r['base'] for r in rr));result={}
            for c in CONTRACTS:
                cc=[r for r in rr if r['contract']==c];result[c]={k:stat([r[k] for r in cc]) for k in ['clearance','objective','test_cost','retained','formation_change','net_value','tests']}
            for old,new in [('D3','D3+'),('D1','D4'),('D1','D3'),('D4','D5'),('D3+','D5'),('D3','D5')]:result[old+'_to_'+new]=stat([by[b,old]['objective']-by[b,new]['objective'] for b in bs])
            result['unavailable_interior_trucks']=sum(t is None for r in rr if r['contract']=='D3+' for t in r['te']);result['truck_count']=3*len(bs);summary['groups'][group][ch]=result
    keys=sorted(set((r['base'],r['channel'],r['sample'],r['duration']) for r in costs));matched=[];regions=[];freq=collections.Counter();ties=collections.Counter()
    for key in keys:
        rr=[r for r in costs if (r['base'],r['channel'],r['sample'],r['duration'])==key];fixed=[r for r in rr if r['contract']!='joint'];best=min(fixed,key=lambda r:r['objective']);ex=next(r for r in rr if r['contract']=='joint' and r['method']=='exact');my=next(r for r in rr if r['contract']=='joint' and r['method']=='myopic');assert ex['objective']<=best['objective']+1e-7 and ex['objective']<=my['objective']+1e-7
        label='stop' if best['tests']==0 else best['contract'];freq[label]+=1;co=[r['contract'] for r in fixed if abs(r['objective']-best['objective'])<1e-7];ties[len(co)]+=1
        regions.append(dict(base=key[0],group=best['group'],channel=key[1],sample=key[2],duration=key[3],preferred=label,cooptimal=co,best_fixed_objective=best['objective'],best_fixed_clearance=best['clearance'],best_fixed_test_cost=best['test_cost'],best_fixed_net_value=best['net_value'],joint_objective=ex['objective'],joint_first_action=ex['first_action']))
        matched.append(dict(base=key[0],channel=key[1],sample=key[2],duration=key[3],joint_over_fixed=best['objective']-ex['objective'],myopic_regret=my['objective']-ex['objective']))
    summary['costs']=dict(cases=len(keys),preferred_canonical=dict(freq),cooptimal_count_distribution=dict(ties),joint_over_fixed=stat([r['joint_over_fixed'] for r in matched]),myopic_regret=stat([r['myopic_regret'] for r in matched]))
    coverage=json.loads((ROOT/'results/coverage.json').read_text())['rows'];cov=[]
    for c,mode in itertools.product(CONTRACTS,['normal','forced_first','oracle_location','no_query']):
        cc=[r for r in coverage if r['contract']==c and r['mode']==mode];bad=[r['omitted'] for r in cc];q=sum(x['queried'] for x in bad);detections=sum(x['queried'] and x['rejected'] for x in bad);ad=[x for x in bad if not x['rejected']]
        cov.append(dict(contract=c,mode=mode,n=len(cc),queried=q,query_fraction=q/len(cc),detected=detections,conditional_detection=None if not q else detections/q,admitted=len(ad),admission_fraction=len(ad)/len(cc),supported_rejected=sum(r['supported']['rejected'] for r in cc),deficits=sum(x['reserve_deficit']>1e-7 for x in ad),overlaps=sum(x['overlap']>1e-7 for x in ad),worst_deficit=max([x['reserve_deficit'] for x in ad],default=None),worst_overlap=max([x['overlap'] for x in ad],default=None)))
    cont=json.loads((ROOT/'results/continuous.json').read_text())['rows'];contby={(r['box'],r['true_shape'],r['noise_scale'],r['contract']):r for r in cont};cm=[]
    for box,sh,scale in sorted(set((r['box'],r['true_shape'],r['noise_scale']) for r in cont)):
        rr={c:contby[box,sh,scale,c] for c in ['D1','D3','D3+','D4','D5']}
        unavailable = rr['D3+']['unavailable_extra']
        identified = rr['D4'].get('same_support_as_D5', False)
        gain_low, gain_high = closure_difference(rr['D3'], rr['D3+'], unavailable)
        gap_low, gap_high = closure_difference(rr['D4'], rr['D5'], identified)
        cm.append(dict(box=box,shape=sh,noise=scale,D3plus_gain_lower=gain_low,D3plus_gain_upper=gain_high,
                       D4_gap_to_D5_lower=gap_low,D4_gap_to_D5_upper=gap_high,unavailable=unavailable,
                       D4_same_support_as_D5=identified))
    summary['continuous']=dict(cases=len(cont),converged=sum(r['converged'] for r in cont),
                               max_width=max(r['certificate_width'] for r in cont),
                               max_ideal_width=max(r['certificate_width'] for r in cont if r['noise_scale']==0),
                               exact_D4_identified=sum(r['D4_same_support_as_D5'] for r in cm),comparisons=cm)
    # Pairwise break-even for the fixed zero-marginal-price policies. Expected
    # count and duration differences are reconstructed from their full leaves.
    thresholds=[]
    for base in design()['bases']:
        if base['rep']!=0:continue
        for ch in ['ideal','bounded']:
            m=Planner(bm.Physics(base),ch);m.prepare_noinfo();vals={}
            for c in ['D3','D3+']:
                m.configure([c],prices=(.002,0,0));v=m.run('exact',trace=True);mass=sum(x['probability']*len(x['path']) for x in v['leaves']);counts=0;duration=0
                from final_model import duration as schedule
                for x in v['leaves']:
                    for i,contract in x['path']:
                        count,t=schedule(base['states'][i],contract);counts+=x['probability']*count;duration+=x['probability']*t
                vals[c]=dict(net=v['objective'],tests=mass,count=counts,duration=duration)
            delta=vals['D3']['net']-vals['D3+']['net'];dm=vals['D3+']['count']-vals['D3']['count'];dt=vals['D3+']['duration']-vals['D3']['duration']
            thresholds.append(dict(group=base['group'],channel=ch,gross_net_gain=delta,extra_expected_measurements=dm,extra_expected_seconds=dt,G0=m.scale,sample_break_even=delta/(m.scale*dm) if dm>1e-10 else None,scope='Fixed reference policies; duration coefficient zero. Global optimized choices are in the full region table.'))
    for name,obj in [('summary',summary),('cost_regions',regions),('cost_matched',matched),('coverage_summary',cov),('break_even',thresholds)]: (ROOT/'results'/(name+'.json')).write_text(json.dumps(obj,indent=2))
    for name,arr in [('cost_regions',regions),('cost_matched',matched),('coverage_summary',cov),('break_even',thresholds),('continuous_bounds',cont)]:
        with (ROOT/'results'/(name+'.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=sorted(set(k for r in arr for k in r)));w.writeheader();w.writerows(arr)
    print(json.dumps(dict(counts=summary['counts'],costs=summary['costs'],continuous={k:v for k,v in summary['continuous'].items() if k!='comparisons'}),indent=2))
if __name__=='__main__':main()
