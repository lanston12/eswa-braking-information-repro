#!/usr/bin/env python3
"""Independent finite-event derivative hierarchy for ramp/first-order pairs.

On each lag/ramp/stop segment, acceleration magnitude is zero, constant,
linear, or a single exponential approach. Difference of their jerks has at
most one zero (or is identically zero), available analytically. These split
intervals of monotone acceleration difference. Bracket its roots, then bracket
relative-velocity roots on the resulting monotonic intervals. All possible
relative-displacement stationary maxima are therefore candidates. Arithmetic
and root finding are floating point; this is not an interval proof.
"""
import json,math,itertools
from dataclasses import replace
from scipy.optimize import brentq
from benchmark import OUT,Truck,closure
from stress_tests import SmoothTruck,smooth_closure


def description(t,mid):
    if mid<t.lag or mid>t.stop:return ('zero',)
    if isinstance(t,SmoothTruck):return ('smooth',t.a,t.tau,t.lag)
    if mid<t.lag+t.ramp:return ('linear',t.j,t.lag)
    return ('constant',t.a)

def brake(desc,t):
    if desc[0]=='zero':return 0.
    if desc[0]=='constant':return desc[1]
    if desc[0]=='linear':return desc[1]*(t-desc[2])
    return desc[1]*(-math.expm1(-(t-desc[3])/desc[2]))

def jerk_at(desc,t):
    if desc[0] in ['zero','constant']:return 0.,0.
    if desc[0]=='linear':return desc[1],0.
    return desc[1]/desc[2]*math.exp(-(t-desc[3])/desc[2]),1/desc[2]

def shape_closure(lead,follow):
    points=[0.,lead.lag,follow.lag,lead.stop,follow.stop]
    for t in [lead,follow]:
        if isinstance(t,Truck):points.append(t.lag+t.ramp)
    points=sorted(set(points));candidates=list(points)
    velocity=lambda t:follow.state(t)[1]-lead.state(t)[1]
    for lo,hi in zip(points,points[1:]):
        if hi-lo<1e-13:continue
        mid=(lo+hi)/2;ld=description(lead,mid);fd=description(follow,mid)
        gl,rl=jerk_at(ld,mid);gf,rf=jerk_at(fd,mid);jerk_points=[lo,hi]
        if gl>0 and gf>0 and abs(rl-rf)>1e-14:
            root=mid+math.log(gl/gf)/(rl-rf)
            if lo<root<hi:jerk_points.append(root)
        jerk_points.sort();acceleration_points=[lo,hi]
        acc=lambda t:brake(ld,t)-brake(fd,t)
        for a,b in zip(jerk_points,jerk_points[1:]):
            fa=acc(a);fb=acc(b)
            if fa*fb<0:acceleration_points.append(brentq(acc,a,b,xtol=1e-13))
            if abs(fa)<1e-12:acceleration_points.append(a)
            if abs(fb)<1e-12:acceleration_points.append(b)
        acceleration_points=sorted(set(acceleration_points))
        candidates.extend(acceleration_points)
        for a,b in zip(acceleration_points,acceleration_points[1:]):
            fa=velocity(a);fb=velocity(b)
            if fa*fb<0:candidates.append(brentq(velocity,a,b,xtol=1e-13))
    vals=[follow.state(t)[0]-lead.state(t)[0] for t in candidates]
    idx=max(range(len(vals)),key=vals.__getitem__)
    return max(0.,vals[idx]),candidates[idx],len(candidates)


def main():
    configs=json.load(open(OUT/'stress_manifest.json'))['configs'];maximum=0.;rampmaximum=0.;witness=None;count=0
    for cfg in configs:
        hs=[[Truck(**h) for h in pair] for pair in cfg['hypotheses']];external=Truck(**cfg['external']);comm=cfg['comm']
        nominal=[external]+[t for pair in hs for t in pair]
        for li,l in enumerate(nominal):
            for fi,f in enumerate(nominal[1:]):
                for lm,fm in itertools.product([0,1],repeat=2):
                    lead=l if lm==0 else SmoothTruck(l);follow=replace(f,lag=f.lag+comm) if fm==0 else SmoothTruck(f,comm)
                    v,at,ncand=shape_closure(lead,follow);old=smooth_closure(lead,follow,768);err=abs(v-old);count+=1
                    if err>maximum:maximum=err;witness=dict(config=cfg['id'],lead=li,follow=fi,lead_mode=lm,follow_mode=fm,event_cost=v,grid_cost=old,at=at,candidates=ncand)
                    if lm==fm==0:rampmaximum=max(rampmaximum,abs(v-closure(lead,follow)[0]))
    result=dict(pair_mode_comparisons=count,max_abs_difference_vs_root_grid_m=maximum,max_ramp_difference_vs_original_exact_event_oracle_m=rampmaximum,worst_case=witness,
        method='Finite lag/ramp/stop events; analytic jerk-difference zeros, bracketed acceleration-difference zeros, bracketed relative-velocity zeros; evaluate all resulting displacement candidates.',scope='Completeness follows from piecewise analytic derivative structure for these two response families only; floating arithmetic, not interval proof. Frozen campaign artifacts unchanged.')
    (OUT/'independent_shape_oracle_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
