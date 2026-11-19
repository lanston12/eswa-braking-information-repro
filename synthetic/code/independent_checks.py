"""Separate two-query enumeration and high-precision continuous checks."""
import itertools,json,hashlib,math,zipfile
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
import mpmath as mp
from final_model import ROOT,design,Planner,bm,CONTRACTS,feature
from continuous import dist


def mp_value(value):
    if isinstance(value, F):
        return mp.mpf(value.numerator) / value.denominator
    return mp.mpf(str(value))


def distance_oracle(a, onset, duration, shape):
    speed = mp.mpf(25)
    if shape == 'ramp':
        ramp_end_speed = speed-a*duration/2
        return speed*onset + speed*duration-a*duration**2/6+ramp_end_speed**2/(2*a)
    stop = mp.findroot(lambda u: a*(u+duration*mp.expm1(-u/duration))-speed,
                       (speed/a, speed/a+duration))
    return speed*(onset+stop) - a*stop**2/2 + a*duration*(stop+duration*mp.expm1(-stop/duration))


def verification_features(state, contract, te):
    a, onset, duration = (F(str(x)) for x in state[:3])
    shape = state[3]
    features = [a]
    if contract in {'D3','D3+'}:
        features.append(onset+F(9,10)*duration if shape=='ramp' else mp_value(onset)+mp.log(10)*mp_value(duration))
    times = ['.35','.65','.95'] if contract=='D4' else ([str(te)] if contract=='D3+' and te is not None else [])
    for time in times:
        elapsed = max(F(0), (F(time)-onset)/duration)
        features.append(a*min(elapsed,1) if shape=='ramp' else
                        F(0) if elapsed==0 else -mp_value(a)*mp.expm1(-mp_value(elapsed)))
    return features


def exact_vertex_oracle(row):
    """Independent high-precision enumeration of all support vertices."""
    truth, contract = row['truth'], row['contract']
    cap, onset, duration = map(mp_value, truth[:3])
    if contract=='D5':
        return distance_oracle(cap,onset,duration,truth[3])
    timing = contract in {'D3','D3+'}
    values = verification_features(truth,contract,row['te'])[1:]
    times = ['.35','.65','.95'] if contract=='D4' else ([str(row['te'])] if contract=='D3+' and row['te'] is not None else [])
    sample_values = values[1:] if timing else values
    tol = mp.mpf('1e-65')
    candidates = []
    for shape in ('ramp','first_order'):
        if shape=='first_order' and any(value==F(str(truth[0])) for value in sample_values):
            continue
        lmin,lmax = map(mp_value,row['prior']['lag'])
        tmin,tmax = map(mp_value,row['prior'][shape])
        inequalities = [(-1,0,-lmin),(1,0,lmax),(0,-1,-tmin),(0,1,tmax)]
        equalities = []
        if timing:
            equalities.append((1,mp.mpf('.9') if shape=='ramp' else mp.log(10),mp_value(values[0])))
        for time,value in zip(map(mp_value,times),sample_values):
            if value==0:
                inequalities.append((-1,0,-time))
            elif value==F(str(truth[0])):
                inequalities.append((1,1,time))
            else:
                normalized = mp_value(value)/cap
                slope = normalized if shape=='ramp' else -mp.log1p(-normalized)
                equalities.append((1,slope,time))
        for (a,b,c),(d,e,f) in itertools.combinations(inequalities+equalities,2):
            determinant = a*e-b*d
            if abs(determinant)<=tol:
                continue
            point = ((c*e-b*f)/determinant,(a*f-c*d)/determinant)
            if all(i*point[0]+j*point[1]<=k+tol for i,j,k in inequalities) and all(abs(i*point[0]+j*point[1]-k)<=tol for i,j,k in equalities):
                candidates.append(distance_oracle(cap,*point,shape))
    return max(candidates)


def check_continuous(rows):
    mp.mp.dps = 80
    leader = mp.mpf(625)/16 + mp.mpf(25)/20 - mp.mpf(1)/300
    exact_count = witness_count = cell_count = 0
    widths = []
    for row in rows:
        for state in [row['truth']]+([row['feasible_maximizer']] if 'feasible_maximizer' in row else []):
            stop = distance_oracle(*map(mp_value,state[:3]),state[3])
            low,high = dist(*state)
            assert mp.mpf(low)<=stop<=mp.mpf(high),(state,low,float(stop),high)
            widths.append(high-low)
        if row['noise_scale']==0:
            optimum = exact_vertex_oracle(row)-leader
            assert mp.mpf(row['closure_lower'])<=optimum<=mp.mpf(row['closure_upper']),row
            exact_count+=1
        elif row['contract']!='D5':
            report = verification_features(row['truth'],row['contract'],row['te'])
            witness = verification_features(row['feasible_maximizer'],row['contract'],row['te'])
            for actual,compatible,feature_report in zip(report,witness,row['feature_reports'],strict=True):
                step=F(str(feature_report['quantum']))
                index=actual//step if isinstance(actual,F) else int(mp.floor(actual/mp_value(step)))
                assert index==feature_report['cell'],row
                low,high=map(F,feature_report['exact_feature_limits'])
                eta=F(str(feature_report['eta']))
                assert (low,high)==(index*step-eta,(index+1)*step+eta),row
                if isinstance(compatible,F):
                    assert low<=compatible<high,row
                else:
                    assert mp_value(low)<=compatible<mp_value(high),row
                cell_count+=1
            assert all(F(str(lo))<=F(str(x))<=F(str(hi)) for x,(lo,hi) in zip(row['feasible_maximizer'][:3],[row['prior']['a'],row['prior']['lag'],row['prior'][row['feasible_maximizer'][3]]],strict=True)),row
            feasible_closure=distance_oracle(*map(mp_value,row['feasible_maximizer'][:3]),row['feasible_maximizer'][3])-leader
            assert mp.mpf(row['closure_lower'])<=feasible_closure,row
            witness_count+=1
    return dict(continuous_distance_points=len(widths),maximum_analytic_distance_enclosure_width=max(widths),
                exact_continuous_vertex_checks=exact_count,quantized_continuous_witness_checks=witness_count,
                quantized_continuous_feature_checks=cell_count)


def compare_previous_enclosures(rows):
    before = json.loads((ROOT/'audit/ideal_observation_before/continuous.json').read_text())['rows']
    paired = list(zip(rows,before,strict=True))
    for current, previous in paired:
        fields = ('box','true_shape','noise_scale','contract','truth','te')
        assert all(current[field]==previous[field] for field in fields)
    exact = [(current,previous) for current,previous in paired if current['noise_scale']==0]
    max_overstatement = max(previous['closure_lower']-current['closure_upper'] for current,previous in exact)
    assert all(current['closure_upper']<=previous['closure_upper']+1e-9 for current,previous in exact)
    return dict(previous_padded_exact_lower_overshoot_m=max(0.0,max_overstatement),
                maximum_bounded_upper_change_m=max(abs(current['closure_upper']-previous['closure_upper'])
                                                   for current,previous in paired if current['noise_scale']>0))


def independent_value(p,m,prices,contracts):
    n=p.n;full=tuple((1<<len(x))-1 for x in p.states);routes=list(itertools.permutations(range(n)))
    def indices(mask):return [k for k in range(8) if mask>>k&1]
    @lru_cache(None)
    def terminal(state):
        ix=[indices(z) for z in state];cost=[]
        for route in routes:
            val=max(p.first[route[0]][ix[route[0]]])+n*2
            val+=sum(max(p.edge[i,j][a,b] for a in ix[i] for b in ix[j]) for i,j in zip(route,route[1:]));val+=.01*m.scale*bm.inversions(route,m.initial)/3;cost.append(val)
        return min(cost)
    reports={}
    for i,c in itertools.product(range(n),contracts):
        o=m.channels[c];vals={}
        for z,mask in enumerate(o.masks[i]):
            bit=sum(1<<k for k,v in enumerate(mask) if v);prob=float(np.dot(p.priors[i],o.likelihood[i][:,z]))
            if prob>1e-13:vals[bit]=vals.get(bit,0)+prob
        reports[i,c]=list(vals.items())
    def update(s,i,z):return tuple(z if k==i else x for k,x in enumerate(s))
    def available(tested):return [(i,c) for i,c in reports if i not in tested and (c!='D3+' or len(contracts)==1 or p.states[i][0]['te'] is not None)]
    @lru_cache(None)
    def one(s,tested):return min([terminal(s)]+[m.fees[i,c]+sum(q*terminal(update(s,i,z)) for z,q in reports[i,c]) for i,c in available(tested)])
    return min([terminal(full)]+[m.fees[i,c]+sum(q*one(update(full,i,z),(i,)) for z,q in reports[i,c]) for i,c in available(())])
def main():
    rows=[];corners=0
    for base in design()['bases']:
        if base['rep']!=0:continue
        p=bm.Physics(base)
        for ch in ['ideal','bounded']:
            m=Planner(p,ch);m.prepare_noinfo();corners+=sum(o.corner_checks for o in m.channels.values())
            for cs,prices in [(['D3+'],(.005,0,0)),(CONTRACTS[:-1],(.002,.01,.01)),(CONTRACTS[:-1],(.002,.04,0))]:
                m.configure(cs,prices=prices);v=independent_value(p,m,prices,cs);w=m.run('exact')['objective'];assert abs(v-w)<1e-8
                rows.append(dict(base=base['id'],channel=ch,contracts=cs,prices=prices,absolute_difference=abs(v-w)))
        print('independent',base['id'],flush=True)
    cont=json.loads((ROOT/'results/continuous.json').read_text())['rows']
    continuous_checks=check_continuous(cont)
    continuous_checks.update(compare_previous_enclosures(cont))
    old=json.loads((ROOT/'audit/FROZEN_R2_MANIFEST.json').read_text());preserved=0
    for rel,v in old['files'].items():
        p=ROOT.parent/'ESWA_R2'/rel
        if p.exists():assert hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256'],rel;preserved+=1
    archived=0
    with zipfile.ZipFile(ROOT/'archive/FROZEN_R2_REPRODUCIBILITY.zip') as z:
        for rel,v in old['files'].items():
            raw=z.read('ESWA_R2/'+rel);assert len(raw)==v['bytes'] and hashlib.sha256(raw).hexdigest()==v['sha256'],rel;archived+=1
    # Uniform ramp continuum, cap=4, t90=0.9, v=25, tau in [0.2,0.8].
    lo,hi=F(1,5),F(4,5);linear=F(25)*(F(1,2)-F(9,10));quad=-F(4,24)
    gain=linear*(lo-(lo+hi)/2)+quad*(lo*lo-(lo*lo+lo*hi+hi*hi)/3)
    assert gain==F(76,25)
    r=dict(status='PASS',independent_two_query_problems=rows,maximum_value_error=max(x['absolute_difference'] for x in rows),new_boundary_corner_checks=corners,**continuous_checks,continuous_all_150_converged=all(x['converged'] for x in cont),maximum_clearance_certificate_width=max(x['certificate_width'] for x in cont),frozen_R2_files_unchanged=preserved,archived_R2_files_verified=archived,uniform_interior_feature_expected_gain_exact=str(gain),uniform_interior_feature_expected_gain_m=float(gain),scope='Separate direct route/two-stage enumeration and high-precision continuous checks. Continuous safety derives from the enclosure proof. All checks were produced in the same assisted research process.')
    (ROOT/'audit/FINAL_CHECKS.json').write_text(json.dumps(r,indent=2));print({k:v for k,v in r.items() if k!='independent_two_query_problems'})
if __name__=='__main__':main()
