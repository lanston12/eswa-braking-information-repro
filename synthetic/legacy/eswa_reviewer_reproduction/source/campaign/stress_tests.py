#!/usr/bin/env python3
"""Explicitly designed challenges, separate from frozen distributional benchmark."""
import sys,json,itertools,math,time
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
from scipy.optimize import brentq
from benchmark import Model,Truck,closure,make_config,run_one,save_results,product_prob,OUT,METHODS


def cap_for_distance(z,v=25.,j=8.,lag=.3):
    return float(brentq(lambda a:Truck(v,a,j,lag).distance-z,.2,50.))

def distance_config(epsilon=.1,delta=.2,distractor=True):
    K=70.;zs=[[K-1,K+3],[K+1,K-3],[K,K]]
    probs=[[1-epsilon,epsilon],[1-epsilon,epsilon],[.5,.5]]
    if distractor:
        zs += [[K+20-delta,K+20+delta],[K+20,K+20]]
        probs += [[.5,.5],[.5,.5]]
    hs=[[asdict(Truck(25.,cap_for_distance(z),8.,.3,f'T{i}_H{h}')) for h,z in enumerate(pair)] for i,pair in enumerate(zs)]
    return dict(id=f'designed_eps{epsilon}_d{delta}_distractor{distractor}',seed=None,family='designed_complementarity',n_uncertain=len(hs),n_total=len(hs)+1,hypotheses=hs,prior=probs,external=asdict(Truck(25.,cap_for_distance(105.),8.,.3,'external')),comm=0.,margin=2.,information='exact independent finite diagnostic labels; duplicate hypotheses represent known trucks')

class SmoothTruck:
    """Alternative first-order brake buildup, tau=a/j matches initial jerk.

    This plant is deliberately outside the certified finite-ramp model class.
    It is a sensitivity challenge, not a validated real-truck plant.
    """
    def __init__(self,t,extra_lag=0.):
        self.v=t.v;self.a=t.a;self.tau=t.a/t.j;self.lag=t.lag+extra_lag
        self.u_stop=brentq(self.raw_velocity,0.,t.v/t.a+self.tau+1.)
        self.stop=self.lag+self.u_stop
        self.distance=self.raw_position(self.u_stop)
    def raw_velocity(self,u):return self.v-self.a*(u-self.tau*(-math.expm1(-u/self.tau)))
    def raw_position(self,u):return self.v*self.lag+self.v*u-self.a*(u*u/2-self.tau*u+self.tau*self.tau*(-math.expm1(-u/self.tau)))
    def state(self,t):
        if t<=self.lag:return self.v*t,self.v
        if t>=self.stop:return self.distance,0.
        u=t-self.lag;return self.raw_position(u),self.raw_velocity(u)

def smooth_closure(lead,follow,resolution=192):
    end=max(lead.stop,follow.stop);grid=np.unique(np.r_[np.linspace(0,end,resolution),lead.lag,follow.lag,lead.stop,follow.stop])
    fun=lambda t:follow.state(float(t))[1]-lead.state(float(t))[1]
    vals=np.array([fun(t) for t in grid]);cand=list(grid)
    for a,b,fa,fb in zip(grid,grid[1:],vals,vals[1:]):
        if fa*fb<0:cand.append(brentq(fun,a,b,xtol=1e-12))
    return max(0.,max(follow.state(t)[0]-lead.state(t)[0] for t in cand))

def actuator_stress(cfg,budget=3):
    m=Model(cfg);n=m.n;k=m.k;actual=np.zeros((n,n,k,k));first=np.zeros((n,k));convergence=0.
    for i in range(n):
        for h in range(k):
            pair=(SmoothTruck(m.external),SmoothTruck(m.hyp[i][h],m.comm))
            first[i,h]=smooth_closure(*pair)
            convergence=max(convergence,abs(first[i,h]-smooth_closure(*pair,resolution=768)))
            for j in range(n):
                if i==j:continue
                for z in range(k):
                    pair=(SmoothTruck(m.hyp[i][h]),SmoothTruck(m.hyp[j][z],m.comm))
                    actual[i,j,h,z]=smooth_closure(*pair)
                    convergence=max(convergence,abs(actual[i,j,h,z]-smooth_closure(*pair,resolution=768)))
    output=[]
    for method in ['no_query','variance','fixed_greedy','myopic','lookahead2','bellman','fullinfo']:
        fail=collision=exceed=worst=0.;marginmean=0.;nominalfail=0.
        for truth in itertools.product(range(k),repeat=n):
            p=product_prob(truth,m.p);state=m.empty
            if method=='fullinfo':state=tuple(h+1 for h in truth)
            elif method!='no_query':
                for rem in range(budget,0,-1):
                    i=m.action(method,state,rem);s=list(state);s[i]=truth[i]+1;state=tuple(s)
            route=m.fixed_route if method=='fixed_greedy' else m.route(state)
            i=route[0];de=[first[i,truth[i]]-m.first[i,state[i]]]
            nominal=[]
            nominal.append(closure(m.external,replace(m.hyp[i][truth[i]],lag=m.hyp[i][truth[i]].lag+m.comm))[0]-m.first[i,state[i]])
            for i,j in zip(route,route[1:]):
                de.append(actual[i,j,truth[i],truth[j]]-m.edge[i,j,state[i],state[j]])
                t=m.hyp[j][truth[j]]
                nominal.append(closure(m.hyp[i][truth[i]],replace(t,lag=t.lag+m.comm))[0]-m.edge[i,j,state[i],state[j]])
            d=max(de);worst=max(worst,d);fail+=p*(d>1e-8);collision+=p*(d>m.margin+1e-8);exceed+=p*max(0,d);marginmean+=p*(m.margin-d);nominalfail+=p*(max(nominal)>1e-8)
        output.append(dict(method=method,budget=budget if method not in ['no_query','fullinfo'] else 0 if method=='no_query' else n,
            probability_violates_2m_certificate=float(fail),probability_negative_clearance=float(collision),expected_max_gap_underestimate_m=float(exceed),worst_gap_underestimate_m=float(worst),expected_min_clearance_m=float(marginmean),in_model_certificate_failure_probability=float(nominalfail)))
    return dict(id=cfg['id'],plant='first-order braking tau=a/j, nominal linear jerk ramp; unchanged cap/lag, exact diagnostic label but unrecognized model family',root_grid_refinement_abs_difference_m=convergence,results=output)


def prior_challenge(configs):
    rows=[]
    for cfg in configs:
        m=Model(cfg)
        from benchmark import prior_stress
        rows.append(dict(id=cfg['id'],prior_stress=prior_stress(m,[2,3])))
    return rows


def correlation_counterexample():
    p={(0,0):.05,(0,1):.45,(1,0):.45,(1,1):.05}
    initial=sum(pr*bits[1] for bits,pr in p.items())
    conditional=sum(pr*bits[1] for bits,pr in p.items() if bits[0]==0)/sum(pr for bits,pr in p.items() if bits[0]==0)
    return dict(joint_probability={str(k):v for k,v in p.items()},reward='min(1,sum revealed theta_i)',initial_item2_marginal=initial,item2_marginal_given_item1_zero=conditional,increase=conditional-initial,full_support=True,
                scope='Logical full-support correlated-outcome counterexample to adaptive diminishing returns for a coverage reward. Not a truck-actuator validation; fixed-order physical embedding requires separate construction.')


def main():
    configs=[distance_config(e,d) for e,d in [(0.01,.02),(.05,.10),(.10,.20),(.20,.40)]]
    configs+=[distance_config(.1,.2,False)]
    results=[]
    for c in configs:
        r,m=run_one(c,budgets=(2,3));results.append(r)
        print(c['id'],{x['method']:round(x['expected_total'],6) for x in r['results'] if x['budget']==2},flush=True)
    save_results(results,configs,'designed')
    # Parameter matrix fixed before this stress run; no rejection by outcome.
    fs=[make_config(853100+q,5,f) for q in range(4) for f in ['heterogeneous','aligned']]
    (OUT/'stress_manifest.json').write_text(json.dumps({'created_before_evaluation':True,'configs':fs,'first_order_tau':'a/j','prior_log_odds_shifts':[-2,-1,1,2]},indent=2))
    ast=[]
    for c in fs:
        r=actuator_stress(c);ast.append(r);print('ACTUATOR',c['id'],r['results'][-2],flush=True)
    (OUT/'actuator_stress.json').write_text(json.dumps(ast,indent=2))
    p=prior_challenge(fs);(OUT/'prior_stress.json').write_text(json.dumps(p,indent=2))
    (OUT/'correlation_counterexample.json').write_text(json.dumps(correlation_counterexample(),indent=2))
    print('DONE',flush=True)
if __name__=='__main__':main()
