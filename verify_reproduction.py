"""Standard-library checks of frozen evidence, not a new vehicle experiment."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from statistics import mean
import csv,hashlib,json,math

ROOT=Path(__file__).resolve().parent
checks=[]
def check(name,condition,detail):
    checks.append(dict(check=name,status='PASS' if condition else 'FAIL',detail=detail))
def rows(rel):
    with (ROOT/rel).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def close(a,b):return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-11)

manifest=json.loads((ROOT/'PROVENANCE.json').read_text())
bad=[d['path'] for d in manifest if hashlib.sha256((ROOT/d['path']).read_bytes()).hexdigest()!=d['sha256']]
check('Frozen-source integrity',not bad,dict(files=len(manifest),mismatches=bad))

def stopping(v,a,lag,tau):return v*lag+v*v/(2*a)+v*tau/2-a*tau*tau/24
leader=stopping(F(25),F(6),F('0.1'),F('0.6'))
def gap(tau):return stopping(F(25),F(4),F('0.9')-F('0.9')*tau,tau)-leader
spread=gap(F('0.2'))-gap(F('0.8'))
check('Exact equal-t90 witness',spread==F('6.1'),dict(spread_m=float(spread),endpoint_closures_m=[float(gap(F('0.2'))),float(gap(F('0.8')))]))
lo=F('.36')/(F('3.6')-(F('2.88')-F('.2')))
hi=F('.36')/(F('3.6')-(F('2.88')+F('.2')))
check('Figure 3 measurement inversion',lo==F(9,23) and hi==F(9,13),dict(posterior_tau_s=[float(lo),float(hi)],prior_reserved_gap_m=float(gap(F('.2'))+2),posterior_reserved_gap_m=float(gap(lo)+2)))

synthetic=rows('synthetic/results/targeted_rows.csv')
check('Finite-support study scope',len(synthetic)==576 and len({r['base'] for r in synthetic})==24,dict(policy_rows=len(synthetic),fleets=len({r['base'] for r in synthetic})))
check('Cost sensitivity scope',len(rows('synthetic/results/targeted_costs.csv'))==2100,dict(rows=len(rows('synthetic/results/targeted_costs.csv'))))

rep=rows('v4/V4_Fig3_policy_replications(1).csv');summary=rows('v4/V4_Fig3_policy_summary(1).csv')
groups=defaultdict(list)
for r in rep:groups[r['scenario'],r['policy']].append(r)
bad=[]
for r in summary:
    g=groups[r['scenario'],r['policy']]
    if len(g)!=int(r['n_replications']) or not close(mean(float(x['heldout_net_gain_m']) for x in g),r['mean_heldout_net_gain_m']):bad.append((r['scenario'],r['policy']))
check('V4 means from individual replications',not bad,dict(groups=len(groups),mismatches=bad))
overlap='Constructed overlapping-response case'
da={r['replication']:float(r['heldout_net_gain_m']) for r in groups[overlap,'Decision-aware D3+']}
ra={r['replication']:float(r['heldout_net_gain_m']) for r in groups[overlap,'Random expected (per-query stop)']}
advantage=mean(da[k]-ra[k] for k in da)
check('V4 primary paired advantage',len(da)==30 and close(advantage,.166287785415103),dict(replications=len(da),advantage_m_equivalent=advantage))
controls=[r for r in rep if r['scenario']!=overlap and r['policy']=='Decision-aware D3+']
check('V4 negative controls retain STOP',len(controls)==60 and all(r['selected_action']=='STOP' and float(r['heldout_net_gain_m'])==0 for r in controls),dict(rows=len(controls)))
forced=[float(r['heldout_net_gain_m']) for r in groups[overlap,'Variance-first (forced)']]
check('V4 adverse forced-query result retained',mean(forced)<0,dict(mean_m_equivalent=mean(forced)))

rep5=rows('v5_3/V5_3_policy_replications.csv');sum5=rows('v5_3/V5_3_policy_convergence.csv');bad=[]
for r in sum5:
    g=[x for x in rep5 if x['scenario']==r['scenario'] and x['method']==r['method']]
    if len(g)!=int(r['n_replications']) or not close(mean(float(x['paired_advantage_mEq']) for x in g),r['paired_advantage_mean_mEq']):bad.append((r['scenario'],r['method']))
check('V5.3 separate resolution experiment',not bad,dict(groups=len(sum5),mismatches=bad))

events=rows('empirical/frozen_inputs/pressure_threshold_holdout.csv');summ=rows('empirical/frozen_inputs/pressure_threshold_summary.csv');bad=[]
for r in summ:
    g=[x for x in events if x['pressure_percent']==r['pressure_percent'] and x['model']==r['model']]
    if len(g)!=int(r['events']) or not close(mean(float(x['first_second_rmse_m_s']) for x in g),r['first_second_rmse_mean_m_s']) or not close(mean(abs(float(x['distance_to_3_error_m'])) for x in g),r['distance_mae_m']):bad.append((r['pressure_percent'],r['model']))
check('Empirical aggregate-to-event reconciliation',not bad,dict(groups=len(summ),mismatches=bad))
g=[r for r in events if r['pressure_percent']=='5'];paired={}
for r in g:paired.setdefault(r['event'],{})[r['model']]=r
# Supplement N compares M1 (delayed step) and M3 (first-order).
speed_gains=[float(v['delayed_step']['first_second_rmse_m_s'])-float(v['first_order']['first_second_rmse_m_s']) for v in paired.values()]
dist_gains=[abs(float(v['delayed_step']['distance_to_3_error_m']))-abs(float(v['first_order']['distance_to_3_error_m'])) for v in paired.values()]
check('Retained empirical mixed outcomes',len(paired)==6 and sum(x>0 for x in speed_gains)==5 and sum(x>0 for x in dist_gains)==3,dict(comparison='M1 delayed step minus M3 first-order',speed_improvements=sum(x>0 for x in speed_gains),distance_improvements=sum(x>0 for x in dist_gains),speed_gains_m_s=speed_gains))
cam=rows('empirical/frozen_inputs/cambridge_extended_fits.csv')
endpoint=[r for r in cam if r['event']=='figure_18b']
check('Cambridge 18b unavailable endpoint retained',len(endpoint)>0 and all(not r['distance_to_3_error_m'] or r['distance_to_3_error_m'].lower()=='nan' for r in endpoint),dict(rows=len(endpoint)))

report=dict(scope='Frozen-result reconciliation plus exact witness; no new measured data or full Monte Carlo rerun',checks=checks,status='PASS' if all(c['status']=='PASS' for c in checks) else 'FAIL')
(ROOT/'verification_report.json').write_text(json.dumps(report,indent=2)+'\n')
for c in checks:print(c['status'],c['check'],json.dumps(c['detail']))
raise SystemExit(0 if report['status']=='PASS' else 1)
