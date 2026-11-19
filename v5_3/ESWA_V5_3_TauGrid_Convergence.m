function ESWA_V5_3_TauGrid_Convergence()
% ESWA_V5_3_TAUGRID_CONVERGENCE
% Focused numerical validation for the ESWA braking-diagnostics paper.
%
% PURPOSE
%   V5 showed that replacing the V4 time-grid/tau-grid evaluator with
%   root-resolved continuous-time endpoint evaluation changes the
%   overlapping-response decision gain by about 0.03 m-equivalent, while
%   the measured time-grid underestimation is only of order 1e-5 m.
%   This script isolates the remaining tau-support discretisation effect.
%
% METHODS
%   - Same three frozen synthetic scenarios as V5.2.
%   - Same 5% feature-error and 20-ms clock-jitter reference condition.
%   - Exact continuous-time maximum of pairwise ramp closure.
%   - Outer-snapped tau posterior supports at N_tau = 21, 51, 101, 201, 401.
%   - Matched planning/test random draws across all N_tau values.
%   - Decision-aware and fair uniform-random-expected policies.
%   - Separate time-grid convergence check.
%
% SCIENTIFIC SCOPE
%   These are synthetic model-based numerical checks, not measured fleet
%   results, not a real sensor validation, and not a continuous-support
%   safety certificate for vehicles outside the declared ramp family.
%
%   The root-resolved method is exact in continuous time up to floating
%   point root solving for the stated ramp equations. Continuous tau support
%   reduction relies on the endpoint-monotonicity condition p >= 0.5,
%   independent product supports, nonnegative ramp onset, and completion
%   of the ramp before the vehicle stops. These conditions are checked.
%
% RUN
%   Extract this package into a new folder and run:
%       ESWA_V5_3_TauGrid_Convergence
%
% MATLAB R2019b+, base MATLAB only.
%
% RUN_MODE:
%   quick    : smoke test only; not manuscript evidence
%   standard : 10 matched repetitions, 100/200 planning/test draws
%   full     : 30 matched repetitions, 250/500 planning/test draws
%
% Outputs:
%   ESWA_V5_3_TAU_GRID_OUTPUT_<MODE>/
%     V5_3_tau_grid_convergence.csv
%     V5_3_policy_convergence.csv
%     V5_3_action_frequency.csv
%     V5_3_time_grid_convergence.csv
%     V5_3_factorised_numerical_error.csv
%     Figure_V5_3_TauGrid_Convergence.pdf/.png/.fig
%     Figure_V5_3_Policy_Convergence.pdf/.png/.fig
%     Figure_V5_3_Action_Stability.pdf/.png/.fig
%     Figure_V5_3_Time_vs_Tau_Error.pdf/.png/.fig
%     V5_3_SCIENTIFIC_AUDIT.txt

clc; close all;
rng(20261116,'twister');

RUN_MODE = 'standard';

switch lower(RUN_MODE)
    case 'quick'
        nRep = 2; nPlan = 30; nTest = 60;
        tauGridCounts = [21 51 101];
    case 'standard'
        nRep = 10; nPlan = 100; nTest = 200;
        tauGridCounts = [21 51 101 201 401];
    case 'full'
        nRep = 30; nPlan = 250; nTest = 500;
        tauGridCounts = [21 51 101 201 401];
    otherwise
        error('Unknown RUN_MODE. Use quick, standard or full.');
end

eta = 0.05;              % assumed bounded b(te) feature error, m/s^2
clockHalfWidth = 0.02;   % assumed clock-jitter half-width, s
queryFee = 0.70;         % normalized clearance-equivalent fee
seedBase = 930000;       % matched V5 seed protocol
dtTime = [0.020 0.010 0.005 0.002 0.001];

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
outDir = fullfile(here,['ESWA_V5_3_TAU_GRID_OUTPUT_' upper(RUN_MODE)]);
if exist(outDir,'dir') ~= 7, mkdir(outDir); end

fprintf('\nESWA V5.3: tau-support discretisation convergence\n');
fprintf('Mode=%s | reps=%d | planning=%d | test=%d\n', ...
    RUN_MODE,nRep,nPlan,nTest);
fprintf('N_tau = %s | eta=%.3f m/s^2 | clock=%.3f s | fee=%.3f m-eq\n', ...
    mat2str(tauGridCounts),eta,clockHalfWidth,queryFee);
fprintf('All fleet results are SYNTHETIC physical-model illustrations.\n\n');

models = makeFrozenScenarios();
scenarioNames = {'Overlap','Fitted proxy','Narrow'};
nScenario = numel(models);
nGrid = numel(tauGridCounts);
nMethod = nGrid + 1;  % method 1 is continuous tau; methods 2.. are snapped grids

% Result arrays: [replication, scenario, method].
daNet = zeros(nRep,nScenario,nMethod);
randomNet = zeros(nRep,nScenario,nMethod);
pairedAdv = zeros(nRep,nScenario,nMethod);
selectedAction = zeros(nRep,nScenario,nMethod); % 0=STOP, 1..3=truck
queryExecuted = false(nRep,nScenario,nMethod);
actionAgreement = false(nRep,nScenario,nMethod);

% Diagnostics of support-grid conservatism, accumulated over independent
% test draws and all candidate trucks.
excessRows = {};
er = 0;
policyRows = {};
pr = 0;
actionRows = {};
ar = 0;

% Keep one exact continuous reference action for each replication/scenario.
for ic = 1:nScenario
    model = models{ic};
    validateModel(model);
    fprintf('Precomputing exact edge matrix for %s...\n',scenarioNames{ic});
    baseEdge = buildRootEdgeMatrix(model,model.tauLo,model.tauHi);
    [baseObjective,baseOrder,baseGaps] = chooseFormation( ...
        model,baseEdge,[],0,0);
    model.baseOrder = baseOrder;
    model.baseObjective = baseObjective;
    model.baseGaps = baseGaps;
    models{ic} = model;

    fprintf('  Base order [%s], base objective %.6f m-equivalent\n', ...
        num2str(baseOrder),baseObjective);

    for ir = 1:nRep
        if ir==1 || mod(ir,2)==0 || ir==nRep
            fprintf('  %s: replication %d/%d\n',scenarioNames{ic},ir,nRep);
            drawnow;
        end

        seedPlan = seedBase + ic*100000 + ir*100;
        seedTest = seedPlan + 1;
        planDraw = makeDraws(nPlan,seedPlan);
        testDraw = makeDraws(nTest,seedTest);

        planGain = zeros(nMethod,3);
        testGain = zeros(nMethod,3);

        % For every truck, compute continuous and all snapped-support
        % outcomes from the SAME observation draws.
        for q = 1:3
            [gExactP,gGridP,~,~,~] = queryBank( ...
                model,baseEdge,baseObjective,baseOrder,q,planDraw, ...
                eta,clockHalfWidth,tauGridCounts);
            [gExactT,gGridT,meanExcess,p95Excess,maxExcess] = queryBank( ...
                model,baseEdge,baseObjective,baseOrder,q,testDraw, ...
                eta,clockHalfWidth,tauGridCounts);

            planGain(1,q) = gExactP;
            testGain(1,q) = gExactT;
            for ig = 1:nGrid
                planGain(ig+1,q) = gGridP(ig);
                testGain(ig+1,q) = gGridT(ig);

                er = er+1;
                excessRows(er,:) = {scenarioNames{ic},ir,q,tauGridCounts(ig), ...
                    meanExcess(ig),p95Excess(ig),maxExcess(ig), ...
                    seedPlan,seedTest};
            end
        end

        % The continuous-support result is the reference action.
        [bestExact,qExact] = max(planGain(1,:));
        if bestExact > queryFee
            exactAction = qExact;
        else
            exactAction = 0;
        end

        for im = 1:nMethod
            pg = planGain(im,:);
            tg = testGain(im,:);
            [bestPlan,qBest] = max(pg);

            % Decision-aware policy: select from planning draws only.
            if bestPlan > queryFee
                thisDA = tg(qBest)-queryFee;
                thisAction = qBest;
                didQuery = true;
            else
                thisDA = 0;
                thisAction = 0;
                didQuery = false;
            end

            % Fair uniform-random expected baseline:
            % each truck is selected with probability 1/3, and each
            % candidate uses its OWN planning value for its stop decision.
            randomThis = mean(double(pg>queryFee).*(tg-queryFee));

            daNet(ir,ic,im) = thisDA;
            randomNet(ir,ic,im) = randomThis;
            pairedAdv(ir,ic,im) = thisDA-randomThis;
            selectedAction(ir,ic,im) = thisAction;
            queryExecuted(ir,ic,im) = didQuery;
            actionAgreement(ir,ic,im) = (thisAction==exactAction);

            if im==1
                methodName = 'Continuous tau / root time';
            else
                methodName = sprintf('N_tau=%d / root time',tauGridCounts(im-1));
            end
            pr = pr+1;
            policyRows(pr,:) = {scenarioNames{ic},ir,methodName, ...
                thisDA,randomThis,thisDA-randomThis,thisAction, ...
                didQuery,exactAction,thisAction==exactAction, ...
                bestPlan,queryFee,seedPlan,seedTest};
            ar = ar+1;
            actionRows(ar,:) = {scenarioNames{ic},ir,methodName, ...
                thisAction,didQuery,exactAction,thisAction==exactAction};
        end
    end
end

% Summaries over replications. These CIs describe Monte Carlo variability
% for the fixed synthetic scenario designs only.
summaryRows = {};
sr = 0;
for ic=1:nScenario
    for im=1:nMethod
        xDA = squeeze(daNet(:,ic,im));
        xR = squeeze(randomNet(:,ic,im));
        xD = squeeze(pairedAdv(:,ic,im));
        xA = squeeze(double(actionAgreement(:,ic,im)));
        xQ = squeeze(double(queryExecuted(:,ic,im)));
        [muD,medD,sdD,loD,hiD,posD,zeroD,negD] = statsCI(xD);
        [muA,medA,sdA,loA,hiA] = statsCI(xDA);
        [muR,medR,sdR,loR,hiR] = statsCI(xR);
        [muAgree,~,~,~,~] = statsCI(xA);
        [muExec,~,~,~,~] = statsCI(xQ);
        if im==1
            methodName='Continuous tau / root time';
            nTau=Inf;
        else
            methodName=sprintf('N_tau=%d / root time',tauGridCounts(im-1));
            nTau=tauGridCounts(im-1);
        end
        sr=sr+1;
        summaryRows(sr,:)={scenarioNames{ic},methodName,nTau,nRep, ...
            muA,sdA,loA,hiA,muR,sdR,loR,hiR, ...
            muD,medD,sdD,loD,hiD,posD,zeroD,negD,muAgree,muExec};
    end
end

policySummary = cell2table(summaryRows,'VariableNames', ...
    {'scenario','method','N_tau','n_replications', ...
     'DA_mean_net_mEq','DA_sd_mEq','DA_t95_lower_mEq','DA_t95_upper_mEq', ...
     'RandomExpected_mean_net_mEq','RandomExpected_sd_mEq', ...
     'RandomExpected_t95_lower_mEq','RandomExpected_t95_upper_mEq', ...
     'paired_advantage_mean_mEq','paired_advantage_median_mEq', ...
     'paired_advantage_sd_mEq','paired_t95_lower_mEq','paired_t95_upper_mEq', ...
     'positive_replication_fraction','zero_replication_fraction', ...
     'negative_replication_fraction','action_agreement_with_continuous', ...
     'query_execution_fraction'});

policyRep = cell2table(policyRows,'VariableNames', ...
    {'scenario','replication','method','DA_net_mEq','RandomExpected_net_mEq', ...
     'paired_advantage_mEq','selected_action','query_executed', ...
     'continuous_reference_action','action_agreement', ...
     'best_planning_gain_mEq','query_fee_mEq','planning_seed','test_seed'});

actionTable = cell2table(actionRows,'VariableNames', ...
    {'scenario','replication','method','selected_action', ...
     'query_executed','continuous_reference_action','action_agreement'});

excessTable = cell2table(excessRows,'VariableNames', ...
    {'scenario','replication','queried_truck','N_tau', ...
     'mean_outer_snap_objective_excess_m','p95_outer_snap_objective_excess_m', ...
     'max_outer_snap_objective_excess_m','planning_seed','test_seed'});

writetable(policySummary,fullfile(outDir,'V5_3_policy_convergence.csv'));
writetable(policyRep,fullfile(outDir,'V5_3_policy_replications.csv'));
writetable(actionTable,fullfile(outDir,'V5_3_action_replications.csv'));

actionSummaryRows={}; asr=0;
for ic=1:nScenario
    for im=1:nMethod
        thisMethod=methodLabel(im,tauGridCounts);
        mask=strcmp(actionTable.scenario,scenarioNames{ic}) & ...
             strcmp(actionTable.method,thisMethod);
        a=actionTable.selected_action(mask);
        asr=asr+1;
        actionSummaryRows(asr,:)={scenarioNames{ic},thisMethod, ...
            numel(a),mean(a==1),mean(a==2),mean(a==3),mean(a==0)};
    end
end
actionSummary=cell2table(actionSummaryRows,'VariableNames', ...
    {'scenario','method','n_replications','truck1_selection_fraction', ...
     'truck2_selection_fraction','truck3_selection_fraction','STOP_fraction'});
writetable(actionSummary,fullfile(outDir,'V5_3_action_frequency.csv'));
writetable(excessTable,fullfile(outDir,'V5_3_tau_grid_convergence.csv'));

%% Independent time-grid convergence and numerical-error decomposition.
timeTable = timeGridConvergence(models,scenarioNames,dtTime);
writetable(timeTable,fullfile(outDir,'V5_3_time_grid_convergence.csv'));

% Summarise tau-grid support conservatism by scenario and grid size.
tauSummaryRows = {};
tr = 0;
for ic=1:nScenario
    for ig=1:nGrid
        mask = strcmp(excessTable.scenario,scenarioNames{ic}) & ...
               excessTable.N_tau==tauGridCounts(ig);
        x = excessTable.mean_outer_snap_objective_excess_m(mask);
        y = excessTable.max_outer_snap_objective_excess_m(mask);
        tr=tr+1;
        tauSummaryRows(tr,:)={scenarioNames{ic},tauGridCounts(ig), ...
            mean(x),simpleQuantile(x,0.95),max(y)};
    end
end
tauSummary = cell2table(tauSummaryRows,'VariableNames', ...
    {'scenario','N_tau','mean_objective_excess_m','p95_mean_excess_m', ...
     'maximum_draw_query_excess_m'});
writetable(tauSummary,fullfile(outDir,'V5_3_tau_grid_summary.csv'));

% Factorised diagnostic: time-grid underestimation is measured separately
% from tau-support outer-snap conservatism. These are distinct quantities.
factorRows = {};
fr=0;
for ic=1:nScenario
    for id=1:numel(dtTime)
        mask=strcmp(timeTable.scenario,scenarioNames{ic}) & ...
             abs(timeTable.dt_s-dtTime(id))<1e-12;
        fr=fr+1;
        factorRows(fr,:)={scenarioNames{ic},'time discretisation', ...
            dtTime(id),NaN,mean(timeTable.underestimation_m(mask)), ...
            max(timeTable.underestimation_m(mask))};
    end
    for ig=1:nGrid
        mask=strcmp(tauSummary.scenario,scenarioNames{ic}) & ...
             tauSummary.N_tau==tauGridCounts(ig);
        fr=fr+1;
        factorRows(fr,:)={scenarioNames{ic},'tau support outer snapping', ...
            NaN,tauGridCounts(ig), ...
            tauSummary.mean_objective_excess_m(mask), ...
            tauSummary.maximum_draw_query_excess_m(mask)};
    end
end
factorTable=cell2table(factorRows,'VariableNames', ...
    {'scenario','error_source','dt_s','N_tau', ...
     'mean_excess_m','maximum_excess_m'});
writetable(factorTable,fullfile(outDir,'V5_3_factorised_numerical_error.csv'));

%% Figures
plotTauGridConvergence(tauSummary,outDir);
plotPolicyConvergence(policySummary,scenarioNames,tauGridCounts,outDir);
plotActionStability(policySummary,scenarioNames,tauGridCounts,outDir);
plotTimeVsTauError(timeTable,tauSummary,scenarioNames,dtTime,tauGridCounts,outDir);
writeAudit(outDir,RUN_MODE,nRep,nPlan,nTest,tauGridCounts,eta, ...
    clockHalfWidth,queryFee,scenarioNames,policySummary,tauSummary,timeTable);

fprintf('\nV5.3 completed. Output folder:\n%s\n',outDir);
fprintf('All outputs are synthetic model-based numerical validation.\n\n');
end

%% ========================== SCENARIOS ==================================
function models=makeFrozenScenarios()
models=cell(3,1);
models{1}=fleetModel([3.90 3.96 4.02],[.75 .73 .71], ...
    [.53 .53 .53],[.78 .78 .78]);
models{2}=fleetModel([3.711807 3.836969 4.457561], ...
    [.881604 .615461 .643483], ...
    [.795372 .530793 .530794],[.979560 .683846 .714981]);
models{3}=fleetModel([3.90 3.96 4.02],[.75 .73 .71], ...
    [.64 .64 .64],[.655 .655 .655]);
end

function m=fleetModel(a,s,lo,hi)
m=struct('v',25,'p',0.90,'a',a(:)','s',s(:)', ...
    'tauLo',lo(:)','tauHi',hi(:)','delta',0.04,'d0',2.0, ...
    'leadA',6.1,'leadTau',0.38,'leadLag',0.10,'beta',0.012);
roundTol=1e-6;
onset=m.s-m.p.*m.tauHi;
if any(onset < -roundTol)
    error('Materially negative ramp onset in frozen scenario.');
end
ix=onset<0;
if any(ix)
    m.tauHi(ix)=(m.s(ix)-1e-9)./m.p;
end
end

function validateModel(m)
if m.p<0.5
    error('Endpoint monotonicity requires p >= 0.5.');
end
if any(m.tauLo>=m.tauHi)
    error('Empty tau support.');
end
if any(m.s-m.p.*m.tauHi < -1e-10)
    error('Negative ramp onset.');
end
if any(m.v./m.a <= m.tauHi/2)
    error('Ramp may not complete before stop; endpoint condition not applicable.');
end
end

%% ========================== DRAWS / OBSERVATIONS =======================
function d=makeDraws(N,seed)
rng(seed,'twister');
d=struct('truth',rand(N,3),'noise',2*rand(N,3)-1, ...
    'clock',2*rand(N,3)-1);
end

function [meanExact,meanGrid,meanExcess,p95Excess,maxExcess]=queryBank( ...
    m,baseEdge,baseObj,baseOrder,q,draw,eta,w,NtauList)
N=size(draw.truth,1);
nG=numel(NtauList);
gExact=zeros(N,1);
gGrid=zeros(N,nG);
excess=zeros(N,nG);
for k=1:N
    lo=m.tauLo; hi=m.tauHi;
    tauTrue=lo(q)+(hi(q)-lo(q))*draw.truth(k,q);
    te=0.5*((m.s(q)-m.p*m.tauLo(q))+m.s(q));
    if phaseEligible(m.s(q),m.p,m.tauLo(q),te,w)
        y=featureValue(m.a(q),m.p,m.s(q),tauTrue,te+w*draw.clock(k,q)) ...
            + eta*draw.noise(k,q);
        [l,h,ok]=posteriorTau(m.a(q),m.p,m.s(q), ...
            m.tauLo(q),m.tauHi(q),te,eta,w,y);
        if ~ok || tauTrue<l-1e-9 || tauTrue>h+1e-9
            error('Posterior support failed to retain a generated truth.');
        end
        lo(q)=l; hi(q)=h;
    end

    objExact=objectiveAfterQuery(m,baseEdge,lo,hi,baseOrder,baseObj,q);
    gExact(k)=baseObj-objExact;

    for ig=1:nG
        [ls,hs]=snapOutward(m.tauLo(q),m.tauHi(q),lo(q),hi(q),NtauList(ig));
        loS=lo; hiS=hi;
        loS(q)=ls; hiS(q)=hs;
        objGrid=objectiveAfterQuery(m,baseEdge,loS,hiS,baseOrder,baseObj,q);
        gGrid(k,ig)=baseObj-objGrid;
        excess(k,ig)=objGrid-objExact;
        if excess(k,ig)<-1e-7
            error('Outer-snapped support produced a smaller robust objective.');
        end
    end
end
meanExact=mean(gExact);
meanGrid=mean(gGrid,1);
meanExcess=mean(excess,1);
p95Excess=zeros(1,nG); maxExcess=zeros(1,nG);
for ig=1:nG
    p95Excess(ig)=simpleQuantile(excess(:,ig),0.95);
    maxExcess(ig)=max(excess(:,ig));
end
end

function [loS,hiS]=snapOutward(tauMin,tauMax,lo,hi,Ntau)
g=linspace(tauMin,tauMax,Ntau);
iLo=find(g<=lo+1e-12,1,'last');
iHi=find(g>=hi-1e-12,1,'first');
if isempty(iLo), iLo=1; end
if isempty(iHi), iHi=numel(g); end
loS=g(iLo); hiS=g(iHi);
if loS>lo+1e-10 || hiS<hi-1e-10
    error('Outer snapping failed to enclose the continuous posterior interval.');
end
end

function yes=phaseEligible(s,p,tauMin,te,w)
maxLag=s-p*tauMin;
yes=(te-w>maxLag+1e-12) && (te+w<s-1e-12);
end

function [lo,hi,ok]=posteriorTau(a,p,s,tauMin,tauMax,te,eta,w,y)
lo=tauMin; hi=tauMax; ok=false;
if ~phaseEligible(s,p,tauMin,te,w), return; end
D=s-te;
denLo=p-(y-eta)/a;
denHi=p-(y+eta)/a;
if denLo<=0, return; end
lo=max(tauMin,(D-w)/denLo);
if denHi>0
    hi=min(tauMax,(D+w)/denHi);
end
if lo>hi+1e-9, return; end
lo=min(hi,max(tauMin,lo));
hi=max(lo,min(tauMax,hi));
ok=true;
end

function y=featureValue(a,p,s,tau,t)
lambda=s-p*tau;
y=a*min(1,max(0,(t-lambda)/tau));
end

%% ========================== FORMATION OBJECTIVE ========================
function E=buildRootEdgeMatrix(m,lo,hi)
E=NaN(4,4);
for lead=0:3
    for follow=1:3
        if lead==follow, continue; end
        E(lead+1,follow+1)=rootEdge(m,lead,follow,lo,hi);
    end
end
end

function obj=objectiveAfterQuery(m,baseEdge,lo,hi,baseOrder,baseObj,q)
E=baseEdge;
for lead=0:3
    for follow=1:3
        if lead==follow, continue; end
        if lead==q || follow==q
            E(lead+1,follow+1)=rootEdge(m,lead,follow,lo,hi);
        end
    end
end
[obj,~,~]=chooseFormation(m,E,baseOrder,m.beta,baseObj);
end

function [best,chosen,gaps]=chooseFormation(m,E,baseOrder,beta,baseObj)
orders=sortrows(perms(1:3));
best=Inf; chosen=[1 2 3]; gaps=zeros(1,3);
for r=1:size(orders,1)
    order=orders(r,:); prev=0; gg=zeros(1,3);
    for k=1:3
        j=order(k);
        gg(k)=m.d0+E(prev+1,j+1);
        prev=j;
    end
    val=sum(gg);
    if ~isempty(baseOrder) && beta>0
        val=val+beta*baseObj*inversionCount(order,baseOrder)/3;
    end
    if val<best-1e-10
        best=val; chosen=order; gaps=gg;
    end
end
end

function c=rootEdge(m,lead,follow,lo,hi)
if lead==0
    aL=m.leadA; tauL=m.leadTau; lagL=m.leadLag;
else
    aL=m.a(lead); tauL=hi(lead); lagL=m.s(lead)-m.p*tauL;
end
aF=m.a(follow); tauF=lo(follow); lagF=m.s(follow)-m.p*tauF;
c=maxClosureRoot(m.v,aL,tauL,lagL,0,aF,tauF,lagF,m.delta);
end

function closure=maxClosureRoot(v,aL,tauL,lagL,cmdL,aF,tauF,lagF,cmdF)
shiftL=lagL+cmdL; shiftF=lagF+cmdF;
stopL=shiftL+v/aL+tauL/2;
stopF=shiftF+v/aF+tauF/2;
tEnd=max(stopL,stopF);
br=[0 shiftL shiftL+tauL stopL shiftF shiftF+tauF stopF tEnd];
br=unique(sort(br(br>=0 & br<=tEnd)));
cand=br(:)';
for ib=1:numel(br)-1
    L=br(ib); R=br(ib+1);
    if R-L<=1e-13, continue; end
    mid=(L+R)/2;
    cF=speedPoly(v,aF,tauF,shiftF,mid);
    cL=speedPoly(v,aL,tauL,shiftL,mid);
    c=cF-cL;
    scale=max(1,max(abs(c)));
    while numel(c)>1 && abs(c(1))<=1e-12*scale
        c(1)=[];
    end
    if numel(c)>1
        rr=realRootsQuadratic(c);
        for j=1:numel(rr)
            if rr(j)>=L-1e-10 && rr(j)<=R+1e-10
                cand(end+1)=min(R,max(L,rr(j))); %#ok<AGROW>
            end
        end
    end
end
xF=rampPosition(cand,v,aF,tauF,lagF,cmdF);
xL=rampPosition(cand,v,aL,tauL,lagL,cmdL);
closure=max(0,max(xF-xL));
end

function c=speedPoly(v,a,tau,shift,t)
u=t-shift; uStop=v/a+tau/2;
if u<0
    c=[0 0 v];
elseif u>=uStop
    c=[0 0 0];
elseif u<tau
    c=[-a/(2*tau),a*shift/tau,v-a*shift^2/(2*tau)];
else
    c=[0 -a v+a*(shift+tau/2)];
end
end

function rr=realRootsQuadratic(c)
c=c(:)'; scale=max(1,max(abs(c))); tol=1e-12*scale;
while numel(c)>1 && abs(c(1))<=tol, c(1)=[]; end
rr=[];
if numel(c)==3
    A=c(1);B=c(2);C=c(3);
    disc=B^2-4*A*C;
    dTol=1e-12*max(1,B^2+abs(4*A*C));
    if disc < -dTol, return; end
    disc=max(0,disc); sd=sqrt(disc);
    if sd<=tol
        rr=-B/(2*A);
    else
        if B>=0,sgn=1;else,sgn=-1;end
        q=-0.5*(B+sgn*sd);
        if abs(q)>tol
            rr=[q/A C/q];
        else
            rr=[(-B+sd)/(2*A) (-B-sd)/(2*A)];
        end
    end
elseif numel(c)==2
    if abs(c(1))>tol, rr=-c(2)/c(1); end
end
rr=rr(isfinite(rr));
end

function x=rampPosition(t,v,a,tau,lag,cmd)
shift=lag+cmd; u=t-shift; uStop=v/a+tau/2;
x=v*t;
r=(u>0)&(u<=tau);
x(r)=v*t(r)-a*u(r).^3/(6*tau);
p=(u>tau)&(u<uStop);
x(p)=v*t(p)-a/2*(u(p)-tau/2).^2-a*tau^2/24;
st=(u>=uStop);
stopPos=v*(shift+uStop)-a/2*(uStop-tau/2)^2-a*tau^2/24;
x(st)=stopPos;
end

function z=inversionCount(order,baseOrder)
rank=zeros(1,3);
for k=1:3, rank(baseOrder(k))=k; end
r=rank(order); z=0;
for i=1:3
    for j=i+1:3
        z=z+(r(i)>r(j));
    end
end
end

%% ========================== TIME GRID CHECK ============================
function T=timeGridConvergence(models,names,dts)
rows={};r=0;
for ic=1:numel(models)
    m=models{ic};
    edges=[0 1;0 2;0 3;1 2;1 3;2 1;2 3;3 1;3 2];
    for ie=1:size(edges,1)
        lead=edges(ie,1);follow=edges(ie,2);
        lo=m.tauLo;hi=m.tauHi;
        if lead==0
            aL=m.leadA;tauL=m.leadTau;lagL=m.leadLag;
        else
            aL=m.a(lead);tauL=hi(lead);lagL=m.s(lead)-m.p*tauL;
        end
        aF=m.a(follow);tauF=lo(follow);lagF=m.s(follow)-m.p*tauF;
        cExact=maxClosureRoot(m.v,aL,tauL,lagL,0,aF,tauF,lagF,m.delta);
        stopL=lagL+vSafe(m.v,aL,tauL);
        stopF=lagF+m.delta+vSafe(m.v,aF,tauF);
        tEnd=max(stopL,stopF);
        for id=1:numel(dts)
            dt=dts(id); tt=0:dt:tEnd;
            if isempty(tt)||tt(end)<tEnd-1e-12, tt=[tt tEnd]; end
            xF=rampPosition(tt,m.v,aF,tauF,lagF,m.delta);
            xL=rampPosition(tt,m.v,aL,tauL,lagL,0);
            cGrid=max(0,max(xF-xL));
            under=max(0,cExact-cGrid);
            r=r+1;
            rows(r,:)={names{ic},lead,follow,dt,cExact,cGrid,under};
        end
    end
end
T=cell2table(rows,'VariableNames',{'scenario','leader','follower', ...
    'dt_s','root_closure_m','grid_closure_m','underestimation_m'});
end

function z=vSafe(v,a,tau)
z=v/a+tau/2;
end

%% ========================== STATISTICS / FIGURES =======================
function name=methodLabel(im,Ntau)
if im==1
    name='Continuous tau / root time';
else
    name=sprintf('N_tau=%d / root time',Ntau(im-1));
end
end

function [mu,med,sd,lo,hi,pos,zer,neg]=statsCI(x)
x=x(:);n=numel(x);
mu=mean(x);med=median(x);sd=std(x,0);
tc=tCritical95(n-1);
half=tc*sd/sqrt(n);
lo=mu-half;hi=mu+half;
tol=1e-10;
pos=mean(x>tol);zer=mean(abs(x)<=tol);neg=mean(x<-tol);
end

function t=tCritical95(df)
if df<=1,t=12.706;
elseif df==2,t=4.303;
elseif df==3,t=3.182;
elseif df==4,t=2.776;
elseif df==5,t=2.571;
elseif df==6,t=2.447;
elseif df==7,t=2.365;
elseif df==8,t=2.306;
elseif df==9,t=2.262;
elseif df==10,t=2.228;
elseif df<=15,t=2.131;
elseif df<=20,t=2.086;
elseif df<=25,t=2.060;
elseif df<=30,t=2.042;
else,t=1.96;
end
end

function q=simpleQuantile(x,p)
x=sort(x(:));n=numel(x);
if n==0,q=NaN;return;end
if n==1,q=x(1);return;end
z=1+(n-1)*p;lo=floor(z);hi=ceil(z);
if lo==hi,q=x(lo);else,q=x(lo)+(z-lo)*(x(hi)-x(lo));end
end

function plotTauGridConvergence(T,outDir)
scenes=unique(T.scenario,'stable');
fig=figure('Color','w','Position',[80 80 1180 460]);
cols=[0.16 0.39 0.59;0.20 0.57 0.60;0.46 0.52 0.59];
subplot(1,2,1);hold on;
for ic=1:numel(scenes)
    mask=strcmp(T.scenario,scenes{ic});
    S=T(mask,:);
    [x,ix]=sort(S.N_tau);
    y=S.mean_objective_excess_m(ix);
    semilogy(x,max(y,1e-10),'-o','LineWidth',2,'MarkerSize',6,'Color',cols(ic,:));
end
xlabel('Number of points in outer-snapped \tau support');
ylabel('Mean excess robust objective (m)');
title('(a) Mean support-grid conservatism');
legend(scenes,'Location','northeast');grid on;box on;

subplot(1,2,2);hold on;
for ic=1:numel(scenes)
    mask=strcmp(T.scenario,scenes{ic});
    S=T(mask,:);
    [x,ix]=sort(S.N_tau);
    y=S.p95_mean_excess_m(ix);
    semilogy(x,max(y,1e-10),'--s','LineWidth',2,'MarkerSize',6,'Color',cols(ic,:));
end
xlabel('Number of points in outer-snapped \tau support');
ylabel('95th percentile of replicate-query mean excess (m)');
title('(b) Upper-tail support-grid conservatism');
legend(scenes,'Location','northeast');grid on;box on;
addSyntheticFooter(fig);
exportFigure(fig,outDir,'Figure_V5_3_TauGrid_Convergence');
end

function plotPolicyConvergence(S,names,Ntau,outDir)
fig=figure('Color','w','Position',[70 70 1250 760]);
for ic=1:3
    subplot(2,3,ic); hold on;
    mask=strcmp(S.scenario,names{ic});
    T=S(mask,:);
    x=Ntau;
    y=zeros(size(x));lo=zeros(size(x));hi=zeros(size(x));
    for j=1:numel(x)
        q=T.N_tau==x(j);
        y(j)=T.paired_advantage_mean_mEq(q);
        lo(j)=T.paired_t95_lower_mEq(q);
        hi(j)=T.paired_t95_upper_mEq(q);
    end
    errorbar(x,y,y-lo,hi-y,'-o','LineWidth',1.8,'CapSize',6);
    qRef=isinf(T.N_tau);
    if any(qRef)
        yline(T.paired_advantage_mean_mEq(qRef),'k--','Continuous reference');
    end
    yline(0,'k:');xlabel('N_\tau');ylabel('DA - random expected (m-eq.)');
    title(names{ic});grid on;box on;
end
for ic=1:3
    subplot(2,3,ic+3);hold on;
    mask=strcmp(S.scenario,names{ic});
    T=S(mask,:);
    x=Ntau;y=zeros(size(x));z=zeros(size(x));
    for j=1:numel(x)
        q=T.N_tau==x(j);
        y(j)=T.DA_mean_net_mEq(q);
        z(j)=T.RandomExpected_mean_net_mEq(q);
    end
    plot(x,y,'-o','LineWidth',1.8);
    plot(x,z,'--s','LineWidth',1.8);
    qRef=isinf(T.N_tau);
    if any(qRef)
        yline(T.DA_mean_net_mEq(qRef),':','DA continuous');
        yline(T.RandomExpected_mean_net_mEq(qRef),'-.','Random continuous');
    end
    xlabel('N_\tau');ylabel('Held-out net gain (m-eq.)');
    title([names{ic} ' policy net gain']);
    legend({'Decision-aware','Random expected'},'Location','best');
    grid on;box on;
end
addSyntheticFooter(fig);
exportFigure(fig,outDir,'Figure_V5_3_Policy_Convergence');
end

function plotActionStability(S,names,Ntau,outDir)
fig=figure('Color','w','Position',[80 80 1120 420]);
for ic=1:3
    subplot(1,3,ic);hold on;
    mask=strcmp(S.scenario,names{ic});
    T=S(mask,:);
    x=Ntau; agree=zeros(size(x));exec=zeros(size(x));
    for j=1:numel(x)
        q=T.N_tau==x(j);
        agree(j)=T.action_agreement_with_continuous(q);
        exec(j)=T.query_execution_fraction(q);
    end
    plot(x,agree,'-o','LineWidth',2);
    plot(x,exec,'--s','LineWidth',2);
    ylim([0 1.05]);xlabel('N_\tau');ylabel('Fraction');
    title(names{ic});legend({'Action agrees with continuous','Query executed'}, ...
        'Location','best');grid on;box on;
end
addSyntheticFooter(fig);
exportFigure(fig,outDir,'Figure_V5_3_Action_Stability');
end

function plotTimeVsTauError(Ttime,Ttau,names,dts,Ntau,outDir)
fig=figure('Color','w','Position',[80 80 1180 440]);
subplot(1,2,1);hold on;
cols=[0.16 0.39 0.59;0.20 0.57 0.60;0.46 0.52 0.59];
for ic=1:3
    mask=strcmp(Ttime.scenario,names{ic});
    S=Ttime(mask,:);
    means=zeros(size(dts));
    for j=1:numel(dts)
        q=abs(S.dt_s-dts(j))<1e-12;
        means(j)=max(S.underestimation_m(q));
    end
    semilogy(dts,max(means,1e-12),'-o','LineWidth',2,'Color',cols(ic,:));
end
xlabel('Time-grid step (s)');ylabel('Maximum time-grid underestimation (m)');
title('(a) Transient-time discretisation');legend(names,'Location','northwest');
grid on;box on;

subplot(1,2,2);hold on;
for ic=1:3
    mask=strcmp(Ttau.scenario,names{ic});
    S=Ttau(mask,:);
    x=S.N_tau;y=S.mean_objective_excess_m;
    [x,ix]=sort(x);y=y(ix);
    semilogy(x,max(y,1e-10),'-s','LineWidth',2,'Color',cols(ic,:));
end
xlabel('N_\tau');ylabel('Mean outer-snap excess objective (m)');
title('(b) Parameter-support discretisation');legend(names,'Location','northeast');
grid on;box on;
addSyntheticFooter(fig);
exportFigure(fig,outDir,'Figure_V5_3_Time_vs_Tau_Error');
end

function addSyntheticFooter(fig)
annotation(fig,'textbox',[0.18 0.005 0.64 0.025], ...
    'String','SYNTHETIC PHYSICAL-MODEL NUMERICAL VALIDATION | NOT MEASURED FLEET DATA', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',8,'Color',[0.42 0.46 0.50]);
end

function exportFigure(fig,outDir,name)
set(fig,'Color','w','PaperPositionMode','auto');
try
    print(fig,fullfile(outDir,[name '.pdf']),'-dpdf','-bestfit');
catch
    print(fig,fullfile(outDir,[name '.pdf']),'-dpdf');
end
print(fig,fullfile(outDir,[name '.png']),'-dpng','-r250');
savefig(fig,fullfile(outDir,[name '.fig']));
end

function writeAudit(outDir,mode,nRep,nPlan,nTest,Ntau,eta,w,fee, ...
    names,S,Ttau,Ttime)
fid=fopen(fullfile(outDir,'V5_3_SCIENTIFIC_AUDIT.txt'),'w');
fprintf(fid,'ESWA V5.3 tau-support convergence audit\n');
fprintf(fid,'=======================================\n\n');
fprintf(fid,'Run mode: %s\n',mode);
fprintf(fid,'Replications: %d; planning draws: %d; test draws: %d\n',nRep,nPlan,nTest);
fprintf(fid,'N_tau values: %s\n',mat2str(Ntau));
fprintf(fid,'Feature error half-width: %.6f m/s^2\n',eta);
fprintf(fid,'Clock-jitter half-width: %.6f s\n',w);
fprintf(fid,'Diagnostic fee: %.6f m-equivalent\n\n',fee);
fprintf(fid,'All three fleet designs are frozen synthetic scenarios.\n');
fprintf(fid,'Monte Carlo intervals quantify only resampling variability within each fixed design.\n');
fprintf(fid,'They are not real-fleet population intervals.\n\n');
fprintf(fid,'Endpoint monotonicity rationale for fixed-cap/fixed-t90 ramp family:\n');
fprintf(fid,'Let q=t-s and cumulative braking B_tau(t)=integral b_tau(u)du.\n');
fprintf(fid,'During the ramp, dB_tau/dtau = a*(q+p*tau)*(p*tau-q)/(2*tau^2).\n');
fprintf(fid,'On the ramp interval q+p*tau >= 0 and q <= (1-p)*tau.\n');
fprintf(fid,'For p>=0.5, p*tau-q >= (2p-1)*tau >= 0.\n');
fprintf(fid,'During the plateau, dB_tau/dtau = a*(p-0.5) >= 0.\n');
fprintf(fid,'Thus cumulative braking is nondecreasing in tau, velocity and displacement\n');
fprintf(fid,'are nonincreasing in tau, and the independent-support worst edge is attained\n');
fprintf(fid,'at tau_leader=upper and tau_follower=lower, under the checked conditions.\n\n');
fprintf(fid,'The root-resolved transient maximum is computed by evaluating phase boundaries\n');
fprintf(fid,'and all real roots of the piecewise quadratic relative-velocity polynomial.\n');
fprintf(fid,'This is a floating-point analytical root calculation, not interval arithmetic.\n\n');
fprintf(fid,'Key interpretation rules:\n');
fprintf(fid,'1. Outer-snapped tau supports are supersets of the exact posterior intervals.\n');
fprintf(fid,'2. Their robust objectives should be no smaller than the exact-support objective.\n');
fprintf(fid,'3. Any excess is parameter-support discretisation conservatism, not temporal error.\n');
fprintf(fid,'4. Time-grid underestimation is reported separately.\n');
fprintf(fid,'5. D3+ measurement error and clock jitter are assumed sensitivity settings.\n');
fprintf(fid,'6. The three-truck fleets are synthetic compositions, not simultaneous tests.\n');
fprintf(fid,'7. No new measured data or field validation are introduced.\n\n');
fprintf(fid,'Mean paired advantage by scenario and N_tau:\n');
for ic=1:numel(names)
    fprintf(fid,'  %s\n',names{ic});
    mask=strcmp(S.scenario,names{ic});
    for k=1:height(S)
        if mask(k)
            fprintf(fid,'    %s: DA-random %.6f, 95%% t CI [%.6f, %.6f], action agreement %.3f\n', ...
                S.method{k},S.paired_advantage_mean_mEq(k), ...
                S.paired_t95_lower_mEq(k),S.paired_t95_upper_mEq(k), ...
                S.action_agreement_with_continuous(k));
        end
    end
end
fclose(fid);
end
