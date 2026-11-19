"""Independent-support extension of the frozen R1 planner. No new optimizer.
Finite-event ramp/first-order dynamics are reused unchanged from the R0 archive.
"""
from pathlib import Path
import sys,itertools,math,time,json
from functools import lru_cache
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'legacy/eswa_reviewer_reproduction'
if not ARCHIVE.exists():ARCHIVE=ROOT.parent/'ESWA_R1_UPGRADE/R0_SNAPSHOT/ESWA_Reproduction/eswa_reviewer_reproduction'
sys.dont_write_bytecode=True;sys.path.insert(0,str(ARCHIVE/'source/campaign'))
from benchmark import Truck,closure
from stress_tests import SmoothTruck
from independent_shape_oracle import shape_closure

def profile(d,delay=0.):
    t=Truck(d['v'],d['a'],d['j'],d['lag']+delay)
    return t if d.get('shape','ramp')=='ramp' else SmoothTruck(t)
def brake(d,t):
    u=max(0.,t-d['lag']);tau=d['a']/d['j']
    return d['a']*(min(u/tau,1.) if d['shape']=='ramp' else -math.expm1(-u/tau))
def features(d,contract):
    a=d['a']
    if contract=='D0':return []
    if contract=='D1':return [a]
    if contract=='D2':return [a,brake(d,.35)]
    if contract=='D3':return [a,d['lag']+(d['a']/d['j'])*(.9 if d['shape']=='ramp' else math.log(10))]
    if contract=='D4':return [a]+[brake(d,t) for t in [.35,.65,.95]]
    if contract=='D5':return [a,d['j'],d['lag'],int(d['shape']=='first_order')]
    raise ValueError(contract)

class Physics:
    def __init__(self,base):
        tic=time.perf_counter();self.base=base;self.n=base['n'];self.states=base['states'];self.priors=[np.array(p) for p in base['priors']]
        self.profiles=[[profile(d) for d in ss] for ss in self.states]
        self.delayed=[[profile(d,base['comm']) for d in ss] for ss in self.states]
        self.external=[profile(dict(base['external'],shape=shape,lag=base['external']['lag']+lag)) for shape in base['external_shapes'] for lag in base['external_hidden_lag']]
        self.first=[];self.edge={}
        for i in range(self.n):
            self.first.append(np.array([max(shape_closure(l,f)[0] for l in self.external) for f in self.delayed[i]]))
            for j in range(self.n):
                if i!=j:self.edge[i,j]=np.array([[shape_closure(l,f)[0] for f in self.delayed[j]] for l in self.profiles[i]])
        self.seconds=time.perf_counter()-tic

class Observation:
    def __init__(self,p,contract,channel):
        tic=time.perf_counter();self.p=p;self.contract=contract;self.channel=channel;self.cells=[];self.weights=[];self.masks=[];self.likelihood=[]
        for i,states in enumerate(p.states):
            fs=np.array([features(d,contract) for d in states]);k=len(states)
            if contract in ['D0','D5'] or channel=='ideal':
                labels=sorted(set(tuple(np.round(f,12)) for f in fs));mat=np.array([[float(tuple(np.round(f,12))==v) for v in labels] for f in fs]);compat=mat>0
                cells=[dict(exact_features=list(v)) for v in labels]
            else:
                dims=fs.shape[1];eta=np.array([.025 if contract=='D3' and j==1 else .05 for j in range(dims)])
                q=np.array([.1 if contract=='D3' and j==1 else .25 for j in range(dims)])
                labels=set()
                for f in fs:
                    bins=[range(math.floor((x-e)/s),math.floor((x+e)/s)+1) for x,e,s in zip(f,eta,q)]
                    labels.update(itertools.product(*bins))
                labels=sorted(labels);mat=np.zeros((k,len(labels)));compat=np.zeros_like(mat,dtype=bool);cells=[]
                for z,v in enumerate(labels):
                    lo=np.array(v)*q;hi=lo+q
                    width=np.maximum(0.,np.minimum(fs+eta,hi)-np.maximum(fs-eta,lo))
                    mat[:,z]=np.prod(width/(2*eta),axis=1)
                    # Quantizer cells [lo,hi), noise bounds closed. Keep even
                    # compatible zero-likelihood endpoint states for safety.
                    compat[:,z]=np.all((fs+eta>=lo-1e-12)&(fs-eta<hi+1e-12),axis=1)
                    cells.append(dict(bins=list(v),lo=lo.tolist(),hi=hi.tolist(),error=eta.tolist(),upper_open=True))
                # Combine reports only when their likelihood rays AND physical
                # consistency supports coincide. Store every original cell.
                groups={}
                for z in range(mat.shape[1]):
                    w=float(mat[:,z].sum())
                    key=(tuple(compat[:,z]),tuple(np.round(mat[:,z]/w,11))) if w>1e-13 else ('zero',z)
                    groups.setdefault(key,[]).append(z)
                indices=list(groups.values());cells=[dict(equivalent_cells=[cells[z] for z in zz]) for zz in indices]
                compat=np.stack([compat[:,zz[0]] for zz in indices],axis=1);mat=np.stack([mat[:,zz].sum(axis=1) for zz in indices],axis=1)
            assert np.max(abs(mat.sum(axis=1)-1))<1e-8
            self.cells.append(cells);self.likelihood.append(mat);self.weights.append(p.priors[i]@mat)
            self.masks.append([np.ones(k,dtype=bool)]+[compat[:,z] for z in range(mat.shape[1])])
        self.seconds=time.perf_counter()-tic
    def children(self,s,i):
        out=[]
        for y,pr in enumerate(self.weights[i]):
            if pr>1e-13:
                t=list(s);t[i]=y;out.append((float(pr),tuple(t)))
        return out

def inversions(r,initial):
    idx={v:k for k,v in enumerate(initial)};a=[idx[x] for x in r]
    return sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))

class Decision:
    METHODS=['no_testing','random','variance','fixed_greedy','fixed_then_reorder','myopic','lookahead2','rollout','bellman','nonadaptive','system_R1']
    def __init__(self,p,o,cfg):
        tic=time.perf_counter();self.p=p;self.o=o;self.cfg=cfg;self.n=p.n;self.empty=(-1,)*self.n;self.budget=cfg['budget'];self.eligible=tuple(cfg['eligible']);self.units=cfg['resources'];self.system_log={}
        self.first=[np.array([max(p.first[i][m]) for m in o.masks[i]]) for i in range(self.n)]
        self.edge={}
        for i,j in p.edge:
            self.edge[i,j]=np.array([[max(p.edge[i,j][np.ix_(mi,mj)].ravel()) if mi.any() and mj.any() else math.inf for mj in o.masks[j]] for mi in o.masks[i]])
        allroutes=np.array(list(itertools.permutations(range(self.n))),dtype=int)
        f,e=self.edges(self.empty);costs=f[allroutes[:,0]]+e[allroutes[:,:-1],allroutes[:,1:]].sum(axis=1)
        self.initial=tuple(allroutes[np.flatnonzero(costs<=min(costs)+1e-9)[0]].tolist());self.scale=float(min(costs)+self.n*p.base['margin'])
        self.routes=np.array([r for r in allroutes if cfg['recourse']=='full' or inversions(r,self.initial)<=1],dtype=int)
        self.penalties=np.array([cfg['change_weight']*self.scale*inversions(r,self.initial)/max(1,self.n*(self.n-1)/2) for r in self.routes])
        self.fees=np.array(self.units)*cfg['cost']*self.scale
        self.variance=[]
        for ps,pi in zip(p.profiles,p.priors):
            ds=np.array([t.distance for t in ps]);self.variance.append(float(pi@(ds-pi@ds)**2))
        self.setup_seconds=time.perf_counter()-tic
    @lru_cache(None)
    def edges(self,s):
        f=np.array([self.first[i][s[i]+1] for i in range(self.n)]);e=np.zeros((self.n,self.n))
        for i,j in self.edge:e[i,j]=self.edge[i,j][s[i]+1,s[j]+1]
        return f,e
    @lru_cache(None)
    def terminal(self,s,fixed=False):
        f,e=self.edges(s)
        if fixed:
            r=self.initial;return float(f[r[0]]+sum(e[i,j] for i,j in zip(r,r[1:]))+self.n*self.p.base['margin']),0.,r
        costs=f[self.routes[:,0]]+e[self.routes[:,:-1],self.routes[:,1:]].sum(axis=1)+self.n*self.p.base['margin'];net=costs+self.penalties
        z=np.flatnonzero(net<=net.min()+1e-9)[0]
        return float(costs[z]),float(self.penalties[z]),tuple(self.routes[z].tolist())
    def stop(self,s,fixed=False):return sum(self.terminal(s,fixed)[:2])
    def available(self,s,b):return [i for i in self.eligible if s[i]<0 and self.units[i]<=b and self.o.contract!='D0']
    @lru_cache(None)
    def lower(self,s):
        # Incoming-edge relaxation of the all-contract information floor.
        val=self.n*self.p.base['margin']
        for j in range(self.n):
            outcomes=[(1.,s[j]+1)] if s[j]>=0 or j not in self.eligible else [(float(q),y+1) for y,q in enumerate(self.o.weights[j]) if q>1e-13]
            for q,zj in outcomes:
                v=self.first[j][zj]
                for i in range(self.n):
                    if i==j:continue
                    zi=[s[i]+1] if s[i]>=0 or i not in self.eligible else [y+1 for y,pr in enumerate(self.o.weights[i]) if pr>1e-13]
                    v=min(v,float(min(self.edge[i,j][zi,zj])))
                val+=q*v
        return float(val)
    @lru_cache(None)
    def value(self,s,b,depth=-1,fixed=False):
        cur=self.stop(s,fixed)
        if b<=0 or depth==0:return cur
        for i in self.available(s,b):cur=min(cur,self.fees[i]+sum(q*self.value(t,b-self.units[i],depth-1 if depth>0 else -1,fixed) for q,t in self.o.children(s,i)))
        return float(cur)
    @lru_cache(None)
    def subset(self,s,b):
        avail=self.available(s,b);best=(self.stop(s),())
        for k in range(1,len(avail)+1):
            for sub in itertools.combinations(avail,k):
                if sum(self.units[i] for i in sub)>b:continue
                def rec(t,rem):
                    if not rem:return self.stop(t)
                    return sum(q*rec(u,rem[1:]) for q,u in self.o.children(t,rem[0]))
                v=sum(self.fees[i] for i in sub)+rec(s,sub)
                if v<best[0]-1e-9:best=(v,sub)
        return best
    @lru_cache(None)
    def action(self,method,s,b):
        avail=self.available(s,b)
        if method=='no_testing' or not avail:return ()
        if method=='system_R1':
            L=self.lower(s)
            if self.stop(s)<=L+min(self.fees[i] for i in avail)+1e-9:self.system_log[(s,b)]='certified_stop';return ()
            U=self.policy_value('myopic',s,b)
            chosen='myopic' if U-L<=.05 else 'bellman' if self.cfg['exact_planning'] else 'lookahead2'
            self.system_log[(s,b)]='bounded_myopic' if chosen=='myopic' else 'exact_small_fleet' if chosen=='bellman' else 'bounded_lookahead'
            return self.action(chosen,s,b)
        fixed=method in ['fixed_greedy','fixed_then_reorder'];S=self.stop(s,fixed)
        depth=-1 if method=='bellman' else 2 if method=='lookahead2' else 1
        vals={}
        for i in avail:
            vals[i]=self.fees[i]+sum(q*(self.policy_value('myopic',t,b-self.units[i]) if method=='rollout' else self.value(t,b-self.units[i],depth-1 if depth>0 else -1,fixed)) for q,t in self.o.children(s,i))
        if method=='random':return tuple(avail) if np.mean(list(vals.values()))<S-1e-9 else ()
        if method=='variance':
            i=max(avail,key=lambda i:(self.variance[i],-i));return (i,) if vals[i]<S-1e-9 else ()
        best=min(vals.values());return (min(i for i in avail if vals[i]<=best+1e-9),) if best<S-1e-9 else ()
    @lru_cache(None)
    def policy_value(self,method,s,b):
        aa=self.action(method,s,b)
        if not aa:return self.stop(s,method=='fixed_greedy')
        return sum(self.fees[i]+sum(q*self.policy_value(method,t,b-self.units[i]) for q,t in self.o.children(s,i)) for i in aa)/len(aa)
    def evaluate(self,method,trace=False):
        sub=self.subset(self.empty,self.budget)[1] if method=='nonadaptive' else None;leaves=[]
        def acts(s,b):return tuple(i for i in sub if s[i]<0)[:1] if sub is not None else self.action(method,s,b)
        @lru_cache(None)
        def rec(s,b):
            aa=acts(s,b)
            if not aa:
                clearance,change,r=self.terminal(s,method=='fixed_greedy');fraction=np.mean([np.mean(self.o.masks[i][s[i]+1]) for i in range(self.n)])
                return np.array([clearance,0.,change,0.,fraction,float(r!=self.initial)])
            return sum(sum(q*rec(t,b-self.units[i]) for q,t in self.o.children(s,i))+np.array([0.,self.fees[i],0.,1.,0.,0.]) for i in aa)/len(aa)
        v=rec(self.empty,self.budget)
        if trace:
            def walk(s,b,pr,cost,path):
                aa=acts(s,b)
                if not aa:
                    g,h,r=self.terminal(s,method=='fixed_greedy');leaves.append(dict(state=s,probability=pr,diagnostic_cost=cost,clearance=g,change_cost=h,route=r,path=path));return
                for i in aa:
                    for q,t in self.o.children(s,i):walk(t,b-self.units[i],pr*q/len(aa),cost+self.fees[i],path+[(i,t[i])])
            walk(self.empty,self.budget,1.,0.,[])
        return v,leaves
    def clear(self,full=True):
        for f in [self.value,self.subset,self.action,self.policy_value,self.lower]:f.cache_clear()
        if full:
            self.edges.cache_clear();self.terminal.cache_clear()
