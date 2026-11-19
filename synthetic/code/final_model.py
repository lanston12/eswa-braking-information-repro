"""One-shot diagnostic menu and exact small-fleet baseline. No new optimizer."""
from pathlib import Path
from functools import lru_cache
import itertools,math,json,time,hashlib
import numpy as np
import base_model as bm
ROOT=Path(__file__).resolve().parents[1]
CONTRACTS=['D1','D2','D3','D3+','D4','D5']
def design():
    p=ROOT/'data/final_design.json';assert hashlib.sha256(p.read_bytes()).hexdigest()==p.with_suffix('.sha256').read_text().strip();d=json.loads(p.read_text())
    extra=ROOT/'data/narrow_response_control.json'
    if extra.exists():
        assert hashlib.sha256(extra.read_bytes()).hexdigest()==extra.with_suffix('.sha256').read_text().strip();d['bases']+=json.loads(extra.read_text())['bases']
    return d
def feature(d,c):
    if c=='D3+':return bm.features(d,'D3')+([] if d.get('te') is None else [bm.brake(d,d['te'])])
    return bm.features(d,c)
def resolution(c,n,scale=1):
    return np.array([.025 if c in ['D3','D3+'] and k==1 else .05 for k in range(n)])*scale,np.array([.1 if c in ['D3','D3+'] and k==1 else .25 for k in range(n)])*scale
class Channel:
    def __init__(self,p,c,channel):
        self.contract=c;self.channel=channel;self.cells=[];self.masks=[];self.weights=[];self.likelihood=[];self.corner_checks=0
        for i,ss in enumerate(p.states):
            ff=np.array([feature(d,c) for d in ss]);k=len(ss)
            if channel=='ideal' or c=='D5':
                labels=sorted(set(tuple(np.round(f,12)) for f in ff));mat=np.array([[float(tuple(np.round(f,12))==v) for v in labels] for f in ff]);mask=mat>0;cells=[dict(exact=list(v)) for v in labels]
            else:
                eta,q=resolution(c,ff.shape[1]);labels=sorted(set(v for f in ff for v in itertools.product(*[range(math.floor((x-e)/s),math.floor((x+e)/s)+1) for x,e,s in zip(f,eta,q)])))
                mat=np.zeros((k,len(labels)));mask=np.zeros_like(mat,dtype=bool);cells=[]
                for z,bins in enumerate(labels):
                    lo=np.array(bins)*q;hi=lo+q
                    mat[:,z]=np.prod(np.maximum(0.,np.minimum(ff+eta,hi)-np.maximum(ff-eta,lo))/(2*eta),axis=1)
                    mask[:,z]=np.all((ff+eta>=lo-1e-12)&(ff-eta<hi+1e-12),axis=1)
                    cells.append(dict(bins=list(bins),lo=lo.tolist(),hi=hi.tolist(),eta=eta.tolist()))
                lookup={tuple(x['bins']):z for z,x in enumerate(cells)}
                for ix,f in enumerate(ff):
                    for signs in itertools.product([-1,1],repeat=len(f)):
                        bins=tuple(np.floor((f+eta*np.array(signs))/q).astype(int));assert mask[ix,lookup[bins]];self.corner_checks+=1
            assert np.max(abs(mat.sum(axis=1)-1))<1e-8
            self.cells.append(cells);self.masks.append(mask.T);self.weights.append(p.priors[i]@mat);self.likelihood.append(mat)
    def actual(self,i,d):
        f=np.array(feature(d,self.contract))
        if self.channel=='ideal' or self.contract=='D5':
            for k,c in enumerate(self.cells[i]):
                if np.allclose(f,c['exact'],atol=1e-11,rtol=0):return k
        else:
            _,q=resolution(self.contract,len(f));bins=np.floor(f/q).astype(int).tolist()
            for k,c in enumerate(self.cells[i]):
                if c['bins']==bins:return k
        return None
def duration(states,c):
    smax=max(feature(d,'D3')[1] for d in states);te=states[0]['te']
    return {'D1':(1,3.),'D2':(2,3.35),'D3':(2,3+smax),'D3+':(2 if te is None else 3,3+smax+(te or 0.)),'D4':(4,4.95),'D5':(0,0.)}[c]
class Planner:
    def __init__(self,p,channel):
        tic=time.perf_counter();self.p=p;self.n=p.n;self.empty=(-1,)*p.n;self.channels={c:Channel(p,c,channel) for c in CONTRACTS};self.local=[];self.outcomes={};self.actual_ids={}
        for i in range(p.n):
            masks=[tuple([True]*len(p.states[i]))];ids={masks[0]:-1}
            for c,o in self.channels.items():
                sums={}
                for y,(m,w) in enumerate(zip(o.masks[i],o.weights[i])):
                    key=tuple(m.tolist())
                    if key not in ids:ids[key]=len(masks)-1;masks.append(key)
                    z=ids[key];self.actual_ids[i,c,y]=z
                    if w>1e-13:sums[z]=sums.get(z,0.)+float(w)
                self.outcomes[i,c]=list(sums.items())
            self.local.append([np.array(x) for x in masks])
        self.first=[np.array([max(p.first[i][m]) for m in self.local[i]]) for i in range(p.n)];self.edge={}
        for i,j in p.edge:self.edge[i,j]=np.array([[np.max(p.edge[i,j][np.ix_(a,b)]) for b in self.local[j]] for a in self.local[i]])
        self.routes=list(itertools.permutations(range(p.n)));costs=[self.route_cost(self.empty,r) for r in self.routes];self.initial=self.routes[int(np.argmin(costs))];self.scale=min(costs);self.penalty=[.01*self.scale*bm.inversions(r,self.initial)/max(1,p.n*(p.n-1)/2) for r in self.routes]
        self.seconds=time.perf_counter()-tic;self.fees={};self.allowed=[]
    def route_cost(self,s,r):return float(self.first[r[0]][s[r[0]]+1]+sum(self.edge[i,j][s[i]+1,s[j]+1] for i,j in zip(r,r[1:]))+self.n*self.p.base['margin'])
    @lru_cache(None)
    def terminal(self,s):
        gs=[self.route_cost(s,r) for r in self.routes];net=np.array(gs)+self.penalty;k=int(np.flatnonzero(net<=min(net)+1e-9)[0]);r=self.routes[k]
        return gs[k],self.penalty[k],r,float(np.mean([self.local[i][s[i]+1].mean() for i in range(self.n)]))
    def configure(self,contracts,prices=None,equal=None):
        self.allowed=list(contracts);self.fees={}
        for i in range(self.n):
            for c in contracts:
                m,t=duration(self.p.states[i],c)
                self.fees[i,c]=self.scale*equal if equal is not None else self.scale*(prices[0]+prices[1]*m+prices[2]*t)
        self.action.cache_clear();self.value.cache_clear();self.metrics.cache_clear()
    def available(self,s):
        return [(i,c) for i in range(self.n) if s[i]<0 for c in self.allowed if c!='D3+' or len(self.allowed)==1 or self.p.states[i][0]['te'] is not None]
    def children(self,s,a):
        i,c=a
        for y,w in self.outcomes[a]:
            t=list(s);t[i]=y
            # A noninformative report still consumes the one permitted query.
            # Its local ID must distinguish queried from unqueried below.
            if y==-1:t[i]=self.noinfo_id[i]
            yield w,tuple(t)
    def prepare_noinfo(self):
        self.noinfo_id=[]
        for i in range(self.n):
            self.noinfo_id.append(len(self.local[i])-1);self.local[i].append(self.local[i][0].copy());self.first[i]=np.append(self.first[i],self.first[i][0])
        for i,j in self.edge:
            x=self.edge[i,j];x=np.concatenate([x,x[0:1,:]],axis=0);x=np.concatenate([x,x[:,0:1]],axis=1);self.edge[i,j]=x
    @lru_cache(None)
    def value(self,s,b,method):
        a=self.action(s,b,method)
        return sum(self.terminal(s)[:2]) if a is None else self.fees[a]+sum(w*self.value(t,b-1,method) for w,t in self.children(s,a))
    @lru_cache(None)
    def action(self,s,b,method):
        if b<=0:return None
        stop=sum(self.terminal(s)[:2]);best=stop;chosen=None
        for a in self.available(s):
            v=self.fees[a]+sum(w*(sum(self.terminal(t)[:2]) if method=='myopic' else self.value(t,b-1,method)) for w,t in self.children(s,a))
            if v<best-1e-9:best=v;chosen=a
        return chosen
    @lru_cache(None)
    def metrics(self,s,b,method):
        a=self.action(s,b,method)
        if a is None:
            g,h,r,f=self.terminal(s);return np.array([g,0.,h,0.,f,float(r!=self.initial)])
        return sum(w*self.metrics(t,b-1,method) for w,t in self.children(s,a))+np.array([0.,self.fees[a],0.,1.,0.,0.])
    def run(self,method,budget=2,trace=False):
        tic=time.perf_counter();v=self.metrics(self.empty,budget,method);out=dict(clearance=float(v[0]),test_cost=float(v[1]),route_cost=float(v[2]),objective=float(sum(v[:3])),tests=float(v[3]),retained=float(v[4]),formation_change=float(v[5]),initial_clearance=self.scale,net_value=self.scale-float(sum(v[:3])),first_action=self.action(self.empty,budget,method),seconds=time.perf_counter()-tic)
        if trace:
            leaves=[]
            def walk(s,b,pr,path):
                a=self.action(s,b,method)
                if a is None:
                    g,h,r,f=self.terminal(s);leaves.append(dict(state=s,probability=pr,route=r,clearance=g,path=path));return
                for w,t in self.children(s,a):walk(t,b-1,pr*w,path+[a])
            walk(self.empty,budget,1.,[]);out['leaves']=leaves
        return out
