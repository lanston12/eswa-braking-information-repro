#!/usr/bin/env python3
"""Synthetic finite-hypothesis diagnosis/formation benchmark.

Classical components: held-karp subset dynamic programming (formation),
Bellman finite-horizon recursion (adaptive information acquisition), greedy and
receding two-step lookahead. No algorithmic novelty is claimed.

The primary benchmark uses exact finite-support robust edge maxima. Earlier
rectangular-hull campaigns are retained as explicitly labeled ablations.
A diagnostic exactly reveals its whole hypothesis, without simulation leakage.
N always counts uncertain trucks; one known external lead is additional.
"""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math, sys, time
from pathlib import Path
from functools import lru_cache
from dataclasses import asdict, replace
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'theory_ordering'))
from envelopes import Truck, closure
OUT=Path(__file__).resolve().parent
COST_TOL=1e-9  # metres; choose lowest truck index / lexicographic subset within tolerance
METHODS=['random','variance','fixed_greedy','fixed_greedy_then_reorder','myopic','lookahead2','rollout','bellman','nonadaptive_opt']


def product_prob(outcome, probs):
    return math.prod(probs[i][h] for i,h in enumerate(outcome))

class Model:
    def __init__(self,cfg):
        self.cfg=cfg; self.n=len(cfg['hypotheses']); self.k=len(cfg['hypotheses'][0]); self.base=self.k+1
        self.hyp=[[Truck(**h) for h in hs] for hs in cfg['hypotheses']]
        self.p=np.array(cfg['prior']); self.external=Truck(**cfg['external']); self.comm=cfg['comm'];self.margin=cfg['margin'];self.robustness=cfg.get('robustness','finite')
        self.powers=self.base**np.arange(self.n)
        self.states=np.array(list(itertools.product(range(self.base),repeat=self.n)),dtype=np.int8)
        # Sorting allows direct mixed-radix indexing, least significant digit = truck 0.
        self.states=self.states[np.argsort(self.states @ self.powers)]
        max_observations=cfg.get('max_observations',3)
        known=np.sum(self.states!=0,axis=1)
        self.states=self.states[(known<=max_observations)|(known==self.n)]
        self.state_lookup=np.full(self.base**self.n,-1,dtype=int)
        self.state_lookup[self.states @ self.powers]=np.arange(len(self.states))
        self.empty=(0,)*self.n
        self.variants=[]
        for hs in self.hyp:
            fast=Truck(hs[0].v,max(t.a for t in hs),max(t.j for t in hs),min(t.lag for t in hs))
            slow=Truck(hs[0].v,min(t.a for t in hs),min(t.j for t in hs),max(t.lag for t in hs))
            self.variants.append([(fast,slow)]+[(t,t) for t in hs])
        tic=time.perf_counter(); self.edge=np.zeros((self.n,self.n,self.base,self.base)); self.first=np.zeros((self.n,self.base))
        for i in range(self.n):
            for u in range(self.base):
                self.first[i,u]=closure(self.external,replace(self.variants[i][u][1],lag=self.variants[i][u][1].lag+self.comm))[0]
                for j in range(self.n):
                    if i==j: continue
                    for v in range(self.base):
                        f=self.variants[j][v][1]
                        self.edge[i,j,u,v]=closure(self.variants[i][u][0],replace(f,lag=f.lag+self.comm))[0]
        if self.robustness=='finite':
            for i in range(self.n):
                for u in range(self.base):
                    li=self.hyp[i] if u==0 else [self.hyp[i][u-1]]
                    self.first[i,u]=max(closure(self.external,replace(f,lag=f.lag+self.comm))[0] for f in li)
                    for j in range(self.n):
                        if i==j:continue
                        for v in range(self.base):
                            fj=self.hyp[j] if v==0 else [self.hyp[j][v-1]]
                            self.edge[i,j,u,v]=max(closure(l,replace(f,lag=f.lag+self.comm))[0] for l in li for f in fj)
        self.edge_seconds=time.perf_counter()-tic
        tic=time.perf_counter(); self.costs,self.route_last,self.pred=self.batch_held_karp(); self.formation_seconds=time.perf_counter()-tic
        self.fixed_route=self.canonical_route(self.empty)
        self.fixed_costs=self.route_costs(self.fixed_route)
        self.variance=[]
        for i,hs in enumerate(self.hyp):
            d=np.array([h.distance for h in hs]);mu=d@self.p[i]
            self.variance.append(float(((d-mu)**2)@self.p[i]))

    def index(self,state):
        idx=int(self.state_lookup[int(np.dot(state,self.powers))])
        if idx<0:raise ValueError('Terminal state not precomputed: '+str(state))
        return idx
    def batch_held_karp(self):
        # Standard Bellman-Held-Karp Hamiltonian PATH DP with fixed external start.
        n=self.n; ns=len(self.states); inf=np.inf
        dp=np.full((1<<n,n,ns),inf); pred=np.full((1<<n,n,ns),-1,dtype=np.int8)
        for j in range(n): dp[1<<j,j]=self.first[j,self.states[:,j]]
        for mask in range(1,1<<n):
            for j in range(n):
                if not mask>>j&1 or mask==(1<<j):continue
                prev=mask^(1<<j)
                for i in range(n):
                    if not prev>>i&1:continue
                    cand=dp[prev,i]+self.edge[i,j,self.states[:,i],self.states[:,j]]
                    better=cand<dp[mask,j]-1e-12
                    dp[mask,j,better]=cand[better];pred[mask,j,better]=i
        last=np.argmin(dp[-1],axis=0); cost=dp[-1,last,np.arange(ns)]
        return cost,last,pred
    def route(self,state):
        idx=self.index(state); mask=(1<<self.n)-1;j=int(self.route_last[idx]); r=[]
        while mask:
            r.append(j); p=int(self.pred[mask,j,idx]);mask^=1<<j;j=p
        return tuple(reversed(r))
    def canonical_route(self,state):
        # Single-state Held-Karp with lexicographic path tie-breaking. This is
        # used for the frozen initial route so fixed-path baselines are stable.
        n=self.n;dp={}
        for j in range(n):dp[(1<<j,j)]=(float(self.first[j,state[j]]),(j,))
        for mask in range(1,1<<n):
            for j in range(n):
                if (mask,j) not in dp:continue
                c,p=dp[(mask,j)]
                for k in range(n):
                    if mask>>k&1:continue
                    key=(mask|1<<k,k);v=c+float(self.edge[j,k,state[j],state[k]]);r=p+(k,)
                    old=dp.get(key)
                    if old is None or v<old[0]-COST_TOL or (abs(v-old[0])<=COST_TOL and r<old[1]):dp[key]=(v,r)
        options=[dp[((1<<n)-1,j)] for j in range(n)];best=min(c for c,r in options)
        return min(r for c,r in options if c<=best+COST_TOL)
    def route_costs(self,r):
        c=self.first[r[0],self.states[:,r[0]]].copy()
        for i,j in zip(r,r[1:]):c+=self.edge[i,j,self.states[:,i],self.states[:,j]]
        return c
    def terminal(self,state,fixed=False):return float((self.fixed_costs if fixed else self.costs)[self.index(state)])
    def children(self,state,i):
        for h in range(self.k):
            s=list(state);s[i]=h+1
            yield h,tuple(s)
    @lru_cache(None)
    def bellman(self,state,budget,fixed=False):
        if budget==0:return self.terminal(state,fixed)
        return min(sum(self.p[i,h]*self.bellman(s,budget-1,fixed) for h,s in self.children(state,i)) for i,u in enumerate(state) if u==0)
    @lru_cache(None)
    def greedy_value(self,state,budget):
        if budget==0:return self.terminal(state)
        i=self.action('myopic',state,budget)
        return sum(self.p[i,h]*self.greedy_value(s,budget-1) for h,s in self.children(state,i))
    @lru_cache(None)
    def action(self,method,state,budget):
        avail=[i for i,u in enumerate(state) if u==0]
        if method=='variance':return max(avail,key=lambda i:(self.variance[i],-i))
        if method=='rollout':
            vals={i:sum(self.p[i,h]*self.greedy_value(s,budget-1) for h,s in self.children(state,i)) for i in avail}
            best=min(vals.values())
            return min(i for i in avail if vals[i]<=best+COST_TOL)
        fixed=(method in ['fixed_greedy','fixed_greedy_then_reorder'])
        depth=budget if method=='bellman' else min(budget,2) if method=='lookahead2' else 1
        vals={i:sum(self.p[i,h]*self.bellman(s,depth-1,fixed) for h,s in self.children(state,i)) for i in avail}
        best=min(vals.values())
        return min(i for i in avail if vals[i]<=best+COST_TOL)
    def nonadaptive_subset(self,budget):
        vals=[]
        for subset in itertools.combinations(range(self.n),budget):
            val=0.
            for hs in itertools.product(range(self.k),repeat=budget):
                s=list(self.empty);p=1.
                for i,h in zip(subset,hs):s[i]=h+1;p*=self.p[i,h]
                val+=p*self.terminal(tuple(s))
            vals.append((val,subset))
        best=min(v for v,s in vals)
        return min(((v,s) for v,s in vals if v<=best+COST_TOL), key=lambda vs:vs[1])
    def evaluate(self,method,budget,true_prior=None,terminal_metric=None):
        # Policy depends only on cfg/model prior and observed state. true_prior is
        # passed solely into the integration weights, never into action selection.
        prior=self.p if true_prior is None else np.array(true_prior)
        subset=self.nonadaptive_subset(budget)[1] if method=='nonadaptive_opt' else None
        @lru_cache(None)
        def rec(s,b):
            if b==0:
                return terminal_metric(s,method=='fixed_greedy') if terminal_metric else self.terminal(s,method=='fixed_greedy')
            avail=[i for i,u in enumerate(s) if u==0]
            actions=avail if method=='random' else [next(i for i in subset if s[i]==0)] if subset is not None else [self.action(method,s,b)]
            return sum(sum(prior[i,h]*rec(t,b-1) for h,t in self.children(s,i)) for i in actions)/len(actions)
        return rec(self.empty,budget)
    def fullinfo(self,prior=None):
        p=self.p if prior is None else prior
        return sum(product_prob(hs,p)*self.terminal(tuple(h+1 for h in hs)) for hs in itertools.product(range(self.k),repeat=self.n))


def make_config(seed,n,family='heterogeneous',k=2):
    rng=np.random.default_rng(seed); hs=[];prior=[]
    for i in range(n):
        if family=='ordered':
            # Deliberate easy control: narrow, fully separated cap supports; tiny
            # boundary/communication effects can leave a nonzero diagnostic gain.
            center=3.0+0.55*i; a=[center,center+0.12]; j=[9.,9.];lag=[.2,.2]
        elif family=='aligned':
            lo=rng.uniform(3.2,5.0);hi=lo+rng.uniform(.5,2.1)
            a=[lo,hi];j=[rng.uniform(2.5,4.0),rng.uniform(5.5,9.)];lag=[rng.uniform(.3,.65),rng.uniform(.08,.23)]
        else:
            a=list(rng.uniform(3.0,7.0,k));j=list(rng.uniform(2.0,10.0,k));lag=list(rng.uniform(.06,.65,k))
        if k!=2 and family!='heterogeneous':raise ValueError('nonbinary only heterogeneous')
        hs.append([asdict(Truck(25.,float(a[h]),float(j[h]),float(lag[h]),f'T{i}_H{h}')) for h in range(k)])
        p=rng.dirichlet(np.full(k,2.0));prior.append(p.tolist())
    return dict(id=f'{family}_n{n}_s{seed}',seed=seed,family=family,n_uncertain=n,n_total=n+1,hypotheses=hs,prior=prior,
                external=asdict(Truck(25.,5.2,6.,.2,'known_external')),comm=.1,margin=2.,robustness='finite',information='exact whole-truck finite hypothesis; independent outcomes')


def audit_small(model,budget=2):
    """Independent permutation and explicit full policy-tree enumeration at N=3."""
    assert model.n==3 and model.k==2
    # No reference to batched Held-Karp cost table in this independent terminal.
    def direct(s):
        best=math.inf
        for r in itertools.permutations(range(3)):
            prev=[model.external];c=0
            for i in r:
                if model.robustness=='finite':
                    followers=model.hyp[i] if s[i]==0 else [model.hyp[i][s[i]-1]]
                    c+=max(closure(l,replace(f,lag=f.lag+model.comm))[0] for l in prev for f in followers)
                    prev=followers
                else:
                    variants=model.variants[i][s[i]];f=variants[1]
                    c+=closure(prev[0],replace(f,lag=f.lag+model.comm))[0];prev=[variants[0]]
            best=min(best,c)
        return best
    errs=[abs(direct(tuple(s))-model.terminal(tuple(s))) for s in model.states]
    # Enumerate first query and a separately selected second query for each of
    # the first query's two observed outcomes: N*(N-1)^K complete policy trees.
    vals=[]
    for first in range(3):
        remain=[i for i in range(3) if i!=first]
        for seconds in itertools.product(remain,repeat=2):
            value=0.
            for h in range(2):
                second=seconds[h]
                for z in range(2):
                    s=[0]*3;s[first]=h+1;s[second]=z+1
                    value+=model.p[first,h]*model.p[second,z]*direct(tuple(s))
            vals.append((value,first,seconds))
    brute=min(vals); dp=model.bellman(model.empty,2)
    return dict(terminal_max_abs_error=max(errs),policy_count=len(vals),brute_policy_value=brute[0],bellman_value=dp,bellman_abs_error=abs(brute[0]-dp),best_first=brute[1],best_seconds=brute[2])


def run_one(cfg,budgets=(2,3),stress=False):
    m=Model(cfg);out=dict(id=cfg['id'],seed=cfg['seed'],family=cfg['family'],n_uncertain=m.n,n_total=m.n+1,
      robustness=m.robustness,policy_cost_tolerance_m=COST_TOL,fixed_route_tie_rule='lexicographically smallest within 1e-9 m Held-Karp local comparisons',edge_seconds=m.edge_seconds,formation_seconds=m.formation_seconds,belief_states=len(m.states),
      no_query_excess=m.terminal(m.empty),no_query_total=m.terminal(m.empty)+m.n*m.margin,
      fullinfo_excess=m.fullinfo(),fixed_route=list(m.fixed_route),results=[])
    for b in budgets:
        if b>m.n:continue
        for method in METHODS:
            # Cold policy caches for comparable online-planning cost; all methods
            # receive the same precomputed terminal costs. Separate formation cost.
            m.bellman.cache_clear();m.action.cache_clear();m.greedy_value.cache_clear();tic=time.perf_counter()
            val=m.evaluate(method,b);sec=time.perf_counter()-tic
            out['results'].append(dict(budget=b,method=method,expected_excess=val,expected_total=val+m.n*m.margin,
                policy_seconds=sec,first_query=None if method in ['random','nonadaptive_opt'] else m.action(method,m.empty,b),
                regret_to_fullinfo=val-out['fullinfo_excess']))
    # Full outcome policy rollouts and sum is an independent evaluation check.
    check=[]
    for b in budgets:
        if b>m.n:continue
        for method in ['myopic','lookahead2','rollout','bellman']:
            v=0
            for hs in itertools.product(range(m.k),repeat=m.n):
                s=m.empty
                for rem in range(b,0,-1):
                    i=m.action(method,s,rem);t=list(s);t[i]=hs[i]+1;s=tuple(t)
                v+=product_prob(hs,m.p)*m.terminal(s)
            target=next(r['expected_excess'] for r in out['results'] if r['budget']==b and r['method']==method)
            check.append(abs(v-target))
    out['rollout_max_abs_error']=max(check,default=0.)
    if stress:
        out['prior_stress']=prior_stress(m,budgets)
    return out,m


def prior_stress(m,budgets):
    if m.k!=2:return []
    output=[]
    for shift in [-2.,-1.,1.,2.]:
        # Shift true log-odds relative to the assumed planning prior. Support stays
        # unchanged, hence the robust certificates remain valid in this model.
        p=m.p[:,1]; q=1/(1+np.exp(-(np.log(p/(1-p))+shift)));true=np.stack([1-q,q],axis=1)
        oracle_cfg=dict(m.cfg);oracle_cfg['prior']=true.tolist();oracle=Model(oracle_cfg)
        for b in budgets:
            if b>m.n:continue
            ideal=oracle.bellman(oracle.empty,b)
            for method in ['variance','myopic','lookahead2','rollout','bellman','random']:
                v=m.evaluate(method,b,true_prior=true)
                output.append(dict(logodds_shift=shift,budget=b,method=method,true_expected_excess=v,true_prior_bellman=ideal,regret_to_true_prior_bellman=v-ideal,support_coverage=True))
    return output


def save_results(results,configs,stem):
    (OUT/f'{stem}_configs.json').write_text(json.dumps(configs,indent=2))
    (OUT/f'{stem}_results.json').write_text(json.dumps(results,indent=2))
    flat=[]
    for r in results:
        for x in r['results']:
            flat.append({k:r[k] for k in ['id','family','n_uncertain','n_total','no_query_excess','no_query_total','fullinfo_excess','formation_seconds']}|x)
    if flat:
        with (OUT/f'{stem}_results.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    return summarize(results)


def summarize(results):
    out={}
    for b in sorted({x['budget'] for r in results for x in r['results']}):
        rows=[r for r in results if any(x['budget']==b for x in r['results'])];d={}
        for method in METHODS:
            v=[next(x for x in r['results'] if x['budget']==b and x['method']==method) for r in rows]
            d[method]={'mean_excess':float(np.mean([x['expected_excess'] for x in v])), 'mean_total':float(np.mean([x['expected_total'] for x in v])), 'mean_policy_seconds':float(np.mean([x['policy_seconds'] for x in v])), 'max_policy_seconds':max(x['policy_seconds'] for x in v)}
        gains=[]
        for r in rows:
            a=next(x['expected_excess'] for x in r['results'] if x['budget']==b and x['method']=='myopic')
            z=next(x['expected_excess'] for x in r['results'] if x['budget']==b and x['method']=='bellman')
            gains.append((a-z,r['id']))
        out[b]={'instances':len(rows),'mean_no_query_excess':float(np.mean([r['no_query_excess'] for r in rows])), 'mean_fullinfo_excess':float(np.mean([r['fullinfo_excess'] for r in rows])), 'methods':d,'bellman_strict_wins_vs_myopic':int(sum(g>1e-8 for g,_ in gains)),'mean_bellman_gain_m':float(np.mean([g for g,_ in gains])),'max_bellman_gain':max(gains)}
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['pilot','main','extended','finite'],default='pilot');args=parser.parse_args()
    # Every seed, family and prior is fixed before evaluation. No rejection or
    # instance selection based on algorithm outcomes. Pilot retained separately.
    if args.stage=='pilot':
        configs=[make_config(20261019+q,4,f) for q,f in enumerate(['heterogeneous','aligned','ordered'])]
        audit=audit_small(Model(make_config(104729,3)))
        (OUT/'independent_audit.json').write_text(json.dumps(audit,indent=2));print('AUDIT',json.dumps(audit),flush=True)
    elif args.stage in ['extended','finite']:
        configs=[make_config(742000+100*n+q,n,f) for n in range(4,11) for f in ['heterogeneous','aligned','ordered'] for q in range(12 if n<=8 else 6)]
    else:
        configs=[make_config(631000+100*n+q,n,f) for n in [4,5,6,7,8] for f in ['heterogeneous','aligned','ordered'] for q in range(4)]
        # seed collision across families is intentional paired construction, each
        # configuration is distinct and appears in the manifest.
    if args.stage!='finite':
        for cfg in configs:cfg['robustness']='box'
    manifest=OUT/f'{args.stage}_manifest.json';manifest.write_text(json.dumps({'created_before_evaluation':True,'stage':args.stage,'configuration_count':len(configs),'configs':configs},indent=2))
    results=[]
    for i,cfg in enumerate(configs):
        r,m=run_one(cfg,stress=(args.stage=='pilot'));results.append(r)
        summary=save_results(results,configs,args.stage)
        print(f'PROGRESS {i+1}/{len(configs)} {cfg["id"]} formation_s={r["formation_seconds"]:.4f} '+json.dumps({b:{k:round(next(x['expected_excess'] for x in r['results'] if x['budget']==b and x['method']==k),5) for k in ['myopic','lookahead2','rollout','bellman']} for b in [2,3]}),flush=True)
    (OUT/f'{args.stage}_summary.json').write_text(json.dumps(summary,indent=2));print('SUMMARY',json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
