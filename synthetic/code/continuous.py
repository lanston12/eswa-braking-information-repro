"""Continuous support certificates for exact and quantized observations.

The external leader completes its ramp by every possible follower onset and
has a greater cap. Hence full-transient closure equals stopping-distance
difference. mpmath interval endpoints certify all pruning and reported bounds;
feasible points only supply lower bounds, never safety certificates.
"""
import math,json,heapq,itertools,time
from fractions import Fraction
import numpy as np
from mpmath import iv
from final_model import ROOT,design
from continuous_exact import exact_report, exact_extrema, cell_index, interval as report_interval

iv.dps=50
def interval(x):return iv.mpf(str(float(x)))
def ends(x):return math.nextafter(float(x.a),-math.inf),math.nextafter(float(x.b),math.inf)
def dist(a,l,t,shape):
    return ends(stopping_distance(*map(interval, (a,l,t)), shape))


def stopping_distance(a, lag, duration, shape):
    v = iv.mpf(25)
    if shape == 'ramp':
        return v*lag + v*v/(2*a) + v*duration/2 - a*duration**2/24
    # Stop time u lies in [v/a,v/a+tau]. With E=exp(-u/tau),
    # D=v*l+v^2/(2a)+v*t-a*t^2/2+a*t^2*(E-E^2/2).
    elo = iv.exp(-v/(a*duration)-1)
    ehi = iv.exp(-v/(a*duration))
    base = v*lag + v*v/(2*a) + v*duration - a*duration**2/2
    low = base + a*duration**2*(elo-elo**2/2)
    high = base + a*duration**2*(ehi-ehi**2/2)
    return iv.mpf([low.a, high.b])
LEADER_LO,LEADER_HI=ends(iv.mpf(625)/16+iv.mpf(25)/20-iv.mpf(1)/300)
LEADER=(LEADER_LO+LEADER_HI)/2

def closure_report(low,high):
    # Outward rounding also covers the final binary-float subtraction/addition.
    lower=math.nextafter(low-LEADER_HI,-math.inf);upper=math.nextafter(high-LEADER_LO,math.inf)
    return dict(closure_lower=lower,closure_upper=upper,clearance_upper=math.nextafter(upper+2,math.inf),certificate_width=math.nextafter(upper-lower,math.inf))


def solve_exact(spec, truth, contract):
    tic = time.perf_counter()
    extra_time = spec['te']
    report = exact_report(truth, contract, extra_time)
    cap = report_interval(report.cap)
    shapes = list(exact_extrema(spec, report))
    assert shapes, 'The exact observation must retain its generating state'
    extrema = []
    for point in shapes:
        low, high = ends(stopping_distance(cap, point.lag, point.duration, point.shape))
        record = dict(shape=point.shape, method=point.method, onset=ends(point.lag), duration=ends(point.duration), stopping_bounds=[low, high])
        if point.duration_range is not None:
            record['compatible_duration'] = [ends(point.duration_range[0])[0], ends(point.duration_range[1])[1]]
        extrema.append(record)
    best = max(range(len(shapes)), key=lambda i: extrema[i]['stopping_bounds'][0])
    point = shapes[best]
    bounds = closure_report(extrema[best]['stopping_bounds'][0], max(r['stopping_bounds'][1] for r in extrema))
    return dict(contract=contract, noise_scale=0, te=extra_time, truth=truth,
                unavailable_extra=extra_time is None and contract=='D3+', **bounds,
                converged=bounds['certificate_width']<=.02, nodes=0, interval_rejections=0,
                feature_reports=[dict(enclosure=ends(report_interval(f))) for f in report.feature_values()],
                feasible_maximizer=[float(report.cap), float(point.lag.mid), float(point.duration.mid), point.shape],
                maximizer_representation='Rounded display of the analytic witness in exact_extrema',
                exact_extrema=extrema, same_support_as_D5=len(extrema)==1 and point.method=='unique_truth',
                seconds=time.perf_counter()-tic,
                support_definition='Exact truth features. Each compatible shape is the prior box intersected with affine equations and phase inequalities. Report enclosures are display metadata, not matching tolerances.')


def fiv(a,l,t,shape,c,te):
    a,l,t=map(interval,(a,l,t));out=[a]
    if c in ['D3','D3+']:out.append(l+(iv.mpf('.9') if shape=='ramp' else iv.ln(10))*t)
    times=[te] if c=='D3+' and te is not None else [.35,.65,.95] if c=='D4' else []
    for tm in times:
        u=interval(tm)-l
        if u.b<=0:out.append(iv.mpf(0));continue
        ul,uh=ends(u);u=iv.mpf([max(0,ul),max(0,uh)])
        if shape=='ramp':
            zl,zh=ends(u/t);out.append(a*iv.mpf([min(1,zl),min(1,zh)]))
        else:out.append(a*(1-iv.exp(-u/t)))
    return [ends(x) for x in out]
def feature_bounds(box,shape,c,te):
    lo,hi=box
    slow=fiv(lo[0],hi[1],hi[2],shape,c,te);fast=fiv(hi[0],lo[1],lo[2],shape,c,te)
    out=[(slow[k][0],fast[k][1]) for k in range(len(slow))]
    if c in ['D3','D3+']:
        out[1]=(fiv(lo[0],lo[1],lo[2],shape,'D3',te)[1][0],fiv(hi[0],hi[1],hi[2],shape,'D3',te)[1][1])
    return out
def tighten(box,shape,c,limits):
    lo,hi=[np.array(x,dtype=float) for x in box];lo[0]=max(lo[0],limits[0][0]);hi[0]=min(hi[0],limits[0][1])
    # Outward-rounded linear contraction for the timing measurement.
    if c in ['D3','D3+']:
        k=iv.mpf('.9') if shape=='ramp' else iv.ln(10);s0,s1=map(interval,limits[1])
        lo[1]=max(lo[1],ends(s0-k*interval(hi[2]))[0]);hi[1]=min(hi[1],ends(s1-k*interval(lo[2]))[1])
        lo[2]=max(lo[2],ends((s0-interval(hi[1]))/k)[0]);hi[2]=min(hi[2],ends((s1-interval(lo[1]))/k)[1])
    return lo,hi
def candidates(box,shape,c,te,limits):
    lo,hi=box;mid=(lo+hi)/2;yield mid
    for p in itertools.product(*zip(lo,hi)):yield np.array(p)
    if c in ['D3','D3+']:
        k=.9 if shape=='ramp' else math.log(10)
        for s in [sum(limits[1])/2,limits[1][0]+1e-9,limits[1][1]-1e-9]:
            for a in [lo[0],mid[0],hi[0]]:
                for t in [lo[2],mid[2],hi[2]]:yield np.array([a,s-k*t,t])
                if c=='D3+' and te is not None:
                    for b in [sum(limits[2])/2,limits[2][0]+1e-9,limits[2][1]-1e-9]:
                        if 0<b<a:
                            den=k-b/a if shape=='ramp' else k+math.log1p(-b/a)
                            if den>0:
                                t=(s-te)/den;yield np.array([a,s-k*t,t])
    if c=='D4':
        a=mid[0];ys=[sum(x)/2 for x in limits[1:]]
        for j,k in itertools.combinations(range(3),2):
            if 0<ys[j]<ys[k]<a:
                t=(.3*(k-j))/((ys[k]-ys[j])/a if shape=='ramp' else math.log((a-ys[j])/(a-ys[k])))
                l=[.35,.65,.95][j]-(t*ys[j]/a if shape=='ramp' else -t*math.log1p(-ys[j]/a))
                yield np.array([a,l,t])
def solve(spec,truth,c,scale):
    assert scale>0 and c!='D5'
    tic=time.perf_counter();a,l,t,shape=truth;te=spec['te']
    report=exact_report(truth,c,te)
    limits=[];interior_limits=[];reported=[]
    for k,feature_value in enumerate(report.feature_values()):
        timing=c in ['D3','D3+'] and k==1
        eta=Fraction('.025' if timing else '.05')*scale
        quantum=Fraction('.1' if timing else '.25')*scale
        label=cell_index(feature_value,quantum)
        lower=label*quantum-eta;upper=(label+1)*quantum+eta
        lower_bounds=ends(report_interval(lower));upper_bounds=ends(report_interval(upper))
        # Use the outer bounds to reject boxes and the inner bounds to accept
        # point witnesses. The exact feature support is half-open.
        limits.append((lower_bounds[0],upper_bounds[1]))
        interior_limits.append((lower_bounds[1],upper_bounds[0]))
        reported.append(dict(value=float(report_interval(feature_value).mid), cell=label,
                             admissible_feature=[float(lower),float(upper)],
                             exact_feature_limits=[str(lower),str(upper)], eta=float(eta), quantum=float(quantum)))
    best=dist(a,l,t,shape)[0];bestpoint=truth;heap=[];serial=0;nodes=0;rejected=0
    def add(box,sh):
        nonlocal best,bestpoint,serial,rejected,nodes
        nodes+=1;box=tighten(box,sh,c,limits);lo,hi=box
        if np.any(lo>hi):rejected+=1;return
        ff=feature_bounds(box,sh,c,te)
        if any(u<L or d>H for (d,u),(L,H) in zip(ff,limits)):rejected+=1;return
        upper=dist(lo[0],hi[1],hi[2],sh)[1]
        if upper<=best:return
        for pt in candidates(box,sh,c,te,limits):
            if np.any(pt<lo) or np.any(pt>hi):continue
            ff=fiv(*pt,sh,c,te)
            if all(d>=L and u<=H for (d,u),(L,H) in zip(ff,interior_limits)):
                val=dist(*pt,sh)[0]
                if val>best:best=val;bestpoint=(*map(float,pt),sh)
        if upper>best:serial+=1;heapq.heappush(heap,(-upper,serial,box,sh))
    for sh in ['ramp','first_order']:add((np.array([spec['a'][0],spec['lag'][0],spec[sh][0]]),np.array([spec['a'][1],spec['lag'][1],spec[sh][1]])),sh)
    while heap and -heap[0][0]-best>.02 and nodes<30000:
        neg,_,box,sh=heapq.heappop(heap)
        if -neg<=best:continue
        lo,hi=box;dim=int(np.argmax((hi-lo)*[30,25,25]));mid=(lo[dim]+hi[dim])/2
        if mid==lo[dim] or mid==hi[dim]:heapq.heappush(heap,(neg,serial+1,box,sh));break
        h=hi.copy();h[dim]=mid;ll=lo.copy();ll[dim]=mid;add((lo,h),sh);add((ll,hi),sh)
    upper=max(best,-heap[0][0] if heap else best)
    bounds=closure_report(best,upper)
    return dict(contract=c,noise_scale=scale,te=te,unavailable_extra=te is None and c=='D3+',truth=truth,feature_reports=reported,**bounds,converged=bounds['certificate_width']<=.02,nodes=nodes,interval_rejections=rejected,feasible_maximizer=bestpoint,seconds=time.perf_counter()-tic,support_definition='Union of both original continuous shape boxes intersected with the true half-open feature intervals. Pruning uses closed outer bounds and nontruth witnesses use inner intervals.')
def main():
    cfg=design()['continuous'];rows=[]
    for box in cfg['boxes']:
        te0=box['lag'][1];te1=min(box['lag'][0]+.9*box['ramp'][0],box['lag'][0]+math.log(10)*box['first_order'][0]);box=dict(box,te=(te0+te1)/2 if te0<te1 else None)
        for shape in ['ramp','first_order']:
            truth=(sum(box['a'])/2,sum(box['lag'])/2,sum(box[shape])/2,shape)
            for scale in cfg['noise_scales']:
                for c in cfg['contracts']:
                    if c=='D5':
                        low,high=dist(*truth);r=dict(contract=c,noise_scale=scale,truth=truth,te=box['te'],unavailable_extra=False,**closure_report(low,high),converged=high-low<=.02,nodes=0,seconds=0.)
                    elif scale==0:r=solve_exact(box,truth,c)
                    else:r=solve(box,truth,c,scale)
                    rows.append(dict(box=box['id'],true_shape=shape,prior=box,**r))
                print('continuous',box['id'],shape,scale,flush=True)
    (ROOT/'results/continuous.json').write_text(json.dumps(dict(rows=rows,scope='Continuous prior intervals. Exact reports use affine compatibility and analytic extrema. Quantized reports use outer pruning and independently certified inner witnesses with 50-digit interval arithmetic. D5 is ideal.'),indent=2))
    witness=[]
    for tau in np.linspace(.2,.8,13):
        b=4*(.9-(.9-.81)/tau)
        for eta in [0,.05,.2,.5]:
            tlo=max(.2,4*(.9-.81)/(3.6-b+eta));thi=min(.8,4*(.9-.81)/(3.6-b-eta)) if 3.6-b-eta>0 else .8
            D=lambda z:25*.9+625/8+25*(.5-.9)*z-4*z*z/24
            witness.append(dict(tau=float(tau),sample=b,error=eta,retained_tau=[tlo,thi],D3_closure=D(.2)-(2.5+625/12+7.5-.09),D3plus_closure=D(tlo)-(2.5+625/12+7.5-.09),gain=D(.2)-D(tlo)))
    (ROOT/'results/direct_witness.json').write_text(json.dumps(witness,indent=2))
if __name__=='__main__':main()
