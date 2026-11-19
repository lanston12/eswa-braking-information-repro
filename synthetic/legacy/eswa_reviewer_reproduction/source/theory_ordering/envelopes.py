"""Exact event/candidate braking-envelope oracle (floating arithmetic, not interval proof)."""
from dataclasses import dataclass, replace, asdict
import math
import itertools
import numpy as np

@dataclass(frozen=True)
class Truck:
    v: float
    a: float
    j: float
    lag: float
    name: str = ''

    @property
    def ramp(self):
        return min(self.a/self.j, math.sqrt(2*self.v/self.j))
    @property
    def stop(self):
        r=self.ramp
        return self.lag+r+max(0., self.v-self.j*r*r/2)/self.a
    def state(self, t):
        # return position, velocity, acceleration, jerk at right-hand interior t
        if t < self.lag:
            return self.v*t, self.v, 0., 0.
        u=t-self.lag
        r=self.ramp
        xr=self.v*self.lag + self.v*r-self.j*r**3/6
        vr=max(0.,self.v-self.j*r*r/2)
        if u<r:
            return self.v*self.lag+self.v*u-self.j*u**3/6, self.v-self.j*u*u/2,-self.j*u,-self.j
        z=min(u-r,vr/self.a)
        x=xr+vr*z-self.a*z*z/2
        if t>=self.stop:
            return x,0.,0.,0.
        return x,vr-self.a*z,-self.a,0.
    @property
    def distance(self):
        return self.state(self.stop+1)[0]

def closure(lead,follow):
    """Maximal position closure and argmax, over all t >= 0."""
    bp=sorted(set([0.,lead.lag,follow.lag,lead.lag+lead.ramp,follow.lag+follow.ramp,lead.stop,follow.stop]))
    cand=bp.copy()
    for l,r in zip(bp,bp[1:]):
        if r-l<1e-10: continue
        m=(l+r)/2
        _,vf,af,jf=follow.state(m)
        _,vl,al,jl=lead.state(m)
        c,b,a=vf-vl,af-al,(jf-jl)/2
        if abs(a)<1e-12:
            roots=[] if abs(b)<1e-12 else [-c/b]
        else:
            d=b*b-4*a*c
            roots=[] if d<0 else [(-b-math.sqrt(d))/(2*a),(-b+math.sqrt(d))/(2*a)]
        cand.extend(m+z for z in roots if l-1e-9<=m+z<=r+1e-9)
    vals=[follow.state(t)[0]-lead.state(t)[0] for t in cand]
    k=int(np.argmax(vals))
    return max(0., vals[k]), cand[k]

def matrix(fleet, comm=0.):
    # Pairwise emergency: lead lag is its physical lag; follower gets extra communication lag.
    n=len(fleet)
    C=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i!=j: C[i,j]=closure(fleet[i],replace(fleet[j],lag=fleet[j].lag+comm))[0]
    return C

def score(C,p): return sum(C[i,j] for i,j in zip(p,p[1:]))
def optimum(C):
    return min((score(C,p),p) for p in itertools.permutations(range(len(C))))

def held_karp(C):
    n=len(C); dp={}; pred={}
    for j in range(n): dp[(1<<j,j)]=0.
    for bits in range(1,1<<n):
        for j in range(n):
            key=(bits,j)
            if key not in dp: continue
            for k in range(n):
                if bits>>k&1: continue
                nk=(bits|1<<k,k); cost=dp[key]+C[j,k]
                if cost<dp.get(nk,float('inf')):
                    dp[nk]=cost; pred[nk]=key
    key=min((((1<<n)-1,j) for j in range(n)), key=lambda k: dp[k])
    cost=dp[key]; route=[key[1]]
    while key in pred:
        key=pred[key];route.append(key[1])
    return cost,tuple(route[::-1])
