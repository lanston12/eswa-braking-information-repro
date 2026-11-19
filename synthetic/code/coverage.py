"""Balanced physical support stress: keep sensing and query-location distinct."""
import json
from final_model import ROOT,CONTRACTS,Planner,design,bm
def trace(m,truth,victim,mode):
    s=m.empty;b=2;queried=False;path=[]
    while b:
        if mode=='no_query':break
        a=(victim,m.allowed[0]) if mode in ['forced_first','oracle_location'] and not path else m.action(s,b,'myopic')
        if a is None or (mode=='oracle_location' and path):break
        i,c=a;queried|=i==victim;path.append(i);y=m.channels[c].actual(i,truth[i])
        if y is None:return dict(queried=queried,rejected=True,queries=path,reserve_deficit=None,overlap=None)
        z=m.actual_ids[i,c,y];t=list(s);t[i]=m.noinfo_id[i] if z==-1 else z;s=tuple(t);b-=1
    g,h,r,f=m.terminal(s);deficits=[]
    for k,j in enumerate(r):
        lead=m.p.external[-1] if k==0 else bm.profile(truth[r[k-1]])
        need=bm.shape_closure(lead,bm.profile(truth[j],m.p.base['comm']))[0]
        provided=m.first[j][s[j]+1] if k==0 else m.edge[r[k-1],j][s[r[k-1]]+1,s[j]+1]
        deficits.append(need-provided)
    deficit=max(0.,max(deficits));return dict(queried=queried,rejected=False,queries=path,reserve_deficit=deficit,overlap=max(0.,deficit-m.p.base['margin']),route=r,clearance=g)
def main():
    rows=[];spec=design()
    for base in spec['bases']:
        if base['rep']!=0:continue
        m=Planner(bm.Physics(base),'bounded');m.prepare_noinfo()
        for c in CONTRACTS:
            m.configure([c],equal=spec['equal_cost'])
            for index in spec['coverage']['truth_indices']:
                good=[dict(ss[(index+i)%8]) for i,ss in enumerate(base['states'])]
                for victim in range(base['n']):
                    bad=[dict(x) for x in good];bad[victim]['j']/=2;bad[victim]['lag']+=.4;bad[victim]['shape']='first_order'
                    for mode in spec['coverage']['methods']:
                        supported=trace(m,good,victim,mode);omitted=trace(m,bad,victim,mode)
                        assert not supported['rejected'] and supported['reserve_deficit']<1e-7
                        rows.append(dict(base=base['id'],group=base['group'],truth_index=index,contract=c,victim=victim,mode=mode,supported=supported,omitted=omitted,truth=good,altered=bad))
        print('coverage',base['id'],flush=True)
    n=sum(r['contract']=='D1' and r['mode']=='normal' for r in rows)
    (ROOT/'results/coverage.json').write_text(json.dumps(dict(rows=rows,scope=f'{n} equally weighted synthetic point scenarios per contract/mode; every truck altered once per truth. Four D-r controls are paired with D-w, not additional independent population draws. Detection means empty response support; the oracle knows only which truck to query and is not executable. Zero realized measurement errors, bounded quantized channel.'),indent=2))
if __name__=='__main__':main()
