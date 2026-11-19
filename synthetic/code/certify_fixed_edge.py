"""Exact rational phase polynomials and rational root-isolation intervals."""
import json
import sympy as s
from final_model import ROOT
x=s.Symbol('t');R=s.Rational
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
    v=[u*w for u in a for w in b];return min(v),max(v)
def interval(poly,lr):
    val=(R(0),R(0))
    for c in s.Poly(poly,x).all_coeffs():val=plus(mul(val,lr),(c,c))
    return val
def truck(a,j,lag):
    a,j,lag=map(R,[a,j,lag]);v=R(25);r=lag+a/j;stop=lag+v/a+a/(2*j)
    d=v*lag+v*v/(2*a)+v*a/(2*j)-a**3/(24*j*j)
    return dict(a=a,j=j,lag=lag,r=r,stop=stop,d=d)
def pos(t,mid):
    if mid<=t['lag']:return 25*x
    if mid<t['r']:return 25*x-t['j']*(x-t['lag'])**3/6
    if mid<t['stop']:return 25*x-t['a']*(x-t['lag']-t['a']/(2*t['j']))**2/2-t['a']**3/(24*t['j']**2)
    return t['d']
def close(a,b):
    points=sorted({R(0),*[t[k] for t in [a,b] for k in ['lag','r','stop']]});vals=[(R(0),R(0))]
    for lo,hi in zip(points,points[1:]):
        poly=s.expand(pos(b,(lo+hi)/2)-pos(a,(lo+hi)/2));der=s.diff(poly,x)
        vals.extend([interval(poly,(lo,lo)),interval(poly,(hi,hi))])
        if der!=0:
            for (u,v),multiplicity in s.polys.polytools.intervals(der,eps=R(1,10**25)):
                if v<lo or u>hi:continue
                vals.append(interval(poly,(max(u,lo),min(v,hi))))
    return max(v[0] for v in vals),max(v[1] for v in vals)
def maxi(vals):return max(a for a,b in vals),max(b for a,b in vals)
def half(a):return a[0]/2,a[1]/2
def main():
    # Rounded rational parameters frozen for a symbolic certificate; original
    # unrounded search witness is retained separately.
    pars=[('3.49','8.30','.13'),('3.55','9.10','.39'),('3.69','9.48','.36'),('4.58','2.68','.11')]
    groups=[[truck(a,j,str(R(lag)+d)) for d in [R(0),R('.12')]] for a,j,lag in pars]
    c=[[maxi([close(a,b) for a in groups[i] for b in groups[2+j]]) for j in range(2)] for i in range(2)]
    before=plus(maxi([a for row in c for a in row]),neg(half(plus(maxi(c[0]),maxi(c[1])))))
    after=plus(maxi([c[0][0],c[1][0]]),neg(half(plus(c[0][0],c[1][0]))))
    difference=plus(after,neg(before));assert difference[0]>0
    def ser(v):return {'lower_exact':str(v[0]),'upper_exact':str(v[1]),'lower':float(v[0]),'upper':float(v[1])}
    out=dict(status='PASS',parameters=pars,unobserved_lag_addition=['0','3/25'],edge_cost_intervals=[[ser(a) for a in row] for row in c],initial_marginal=ser(before),conditional_marginal=ser(after),adaptive_violation=ser(difference),root_isolation_width='1e-25',arithmetic='SymPy rational coefficients, real polynomial root isolation, rational interval Horner evaluation; no floating tolerance in proof of positive violation')
    (ROOT/'audit/REEXECUTED_RATIONAL_CERTIFICATE.json').write_text(json.dumps(out,indent=2));print('Certified adaptive increase in',float(difference[0]),float(difference[1]))
if __name__=='__main__':main()
