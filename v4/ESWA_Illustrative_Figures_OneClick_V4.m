function ESWA_Illustrative_Figures_OneClick_V4()
% ESWA_ILLUSTRATIVE_FIGURES_ONECLICK
% Reproducible physics-based SYNTHETIC reference figures for a study of
% diagnostic information sufficiency and robust truck-platoon formation.
%
% Run from MATLAB (R2019b+, base MATLAB, no optional toolboxes):
%     ESWA_Illustrative_Figures_OneClick_V4
%
% SCIENTIFIC USE RULES:
%  * Everything plotted by this file is a controlled synthetic model result.
%    NOTHING in this program is an actual field validation or empirical fit.
%  * The 6.1 m number is an analytical counterexample for a specified ramp
%    response family, not a universal spacing margin.
%  * Measurement error, clock jitter, response supports, sensor prices and
%    synthetic truck fleets are ASSUMPTIONS, not calibrated population data.
%  * A phase-infeasible D3+ measurement is DISABLED, not credited as a query.
%  * Report both positive and zero/negative cases, and never discard cases
%    merely because the results are less attractive.
%  * The fleet experiment evaluates a product-support, per-edge worst-case
%    sum using an outer-snapped parameter grid and a transient time grid.
%    These numerical calculations are NOT a certified continuum safety proof.
%  * Do not publish these illustrative plots as measured/platoon-test results.
%
% Outputs: figures 1-4 in PDF/PNG/FIG, CSV tables and DESIGN_NOTES.txt.
% V4: fair per-query random stopping, exact random-policy expectation,
%     forced variance-first ablation, 10,000-resample paired bootstrap,
%     independent-draw clearance savings by physical edge.
% All V4 outputs remain SYNTHETIC physical-model calculations.

clc; close all;
rng(20261115, 'twister');
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
out = fullfile(here, 'ESWA_V4_FAIR_VALIDATION_OUTPUT');
if exist(out, 'dir') ~= 7, mkdir(out); end

colors.blue  = [0.16 0.39 0.59];
colors.teal  = [0.20 0.57 0.60];
colors.gray  = [0.46 0.52 0.59];
colors.red   = [0.77 0.37 0.31];
colors.green = [0.25 0.57 0.42];
colors.light = [0.92 0.94 0.96];

fprintf('\nESWA: synthetic physics-based illustration. NOT empirical validation.\n');

%% 1. Exact theorem family: known cap and fixed t90, unknown lag/buildup
P.v = 25.0; P.a = 4.0; P.p = 0.9; P.s = 0.9;
P.tauLo = 0.2; P.tauHi = 0.8; P.te = 0.81;
P.d0 = 2.0; P.fee = 0.70;
% Known stronger leader in the exact analytical witness.
P.leadA = 6.0; P.leadTau = 0.6; P.leadLag = 0.1;

assert(P.s-P.p*P.tauHi >= 0, 'Negative onset is not allowed.');
assert(P.v > P.a*P.tauHi/2, 'Ramp must finish before vehicle stops.');

spread = fixedTimeD(P.v,P.a,P.p,P.s,P.tauLo) ...
       - fixedTimeD(P.v,P.a,P.p,P.s,P.tauHi);
assert(abs(spread-6.1) < 1e-10, 'Theorem witness should give 6.1 m.');
fprintf('Analytical D3 information limit: %.8f m (expected 6.1 m).\n',spread);

% Verify the closure difference from the actual complete trajectories.
t = 0:0.003:13;
xLead = rampX(t,P.v,P.leadA,P.leadTau,P.leadLag,0);
xSlow = rampX(t,P.v,P.a,P.tauLo,P.s-P.p*P.tauLo,0);
xFast = rampX(t,P.v,P.a,P.tauHi,P.s-P.p*P.tauHi,0);
Cslow = max([0 max(xSlow-xLead)]);
Cfast = max([0 max(xFast-xLead)]);
assert(abs((Cslow-Cfast)-spread)<0.015, ...
    'Numerical closure does not agree with the analytical information spread.');

etaExample=0.05; clockExample=0.02; truthExample=0.65;
yExample = capFeature(P.v,P.a,P.p,P.s,truthExample,P.te) + 0.01;
[loExample,hiExample,okExample] = posteriorTau(P.a,P.p,P.s, ...
    P.tauLo,P.tauHi,P.te,etaExample,clockExample,yExample);
assert(okExample && loExample<=truthExample && hiExample>=truthExample, ...
    'The held-out constructed truth must remain in the compatible support.');

fig=figure('Color','w','Position',[70 70 1650 650]);
ax1=subplot(1,3,1); set(ax1,'Position',[0.055 0.23 0.255 0.63]); hold on;
tb=0:0.002:1.6;
plot(tb,capFeature(P.v,P.a,P.p,P.s,P.tauLo,tb), ...
    'LineWidth',2.3,'Color',colors.blue);
plot(tb,capFeature(P.v,P.a,P.p,P.s,P.tauHi,tb), ...
    '--','LineWidth',2.3,'Color',colors.teal);
xline(P.s,':','t_{90}','LineWidth',1.2,'Color',colors.gray);
xline(P.te,'-.','t_e','LineWidth',1.2,'Color',colors.red);
xlabel('Time from common command (s)'); ylabel('Brake magnitude (m/s^2)');
title('(a) Same cap and t_{90}');
legend({'Short buildup','Long buildup'},'Location','southeast');
axis([0 1.55 0 4.45]); grid on; box on;

ax2=subplot(1,3,2); set(ax2,'Position',[0.375 0.23 0.255 0.63]); hold on;
plot(t,xSlow-xLead,'LineWidth',2.3,'Color',colors.blue);
plot(t,xFast-xLead,'--','LineWidth',2.3,'Color',colors.teal);
yline(Cslow,':','LineWidth',1.0,'Color',colors.blue);
yline(Cfast,':','LineWidth',1.0,'Color',colors.teal);
xlabel('Time from leader command (s)');
ylabel('Relative travel (m)');
title('(b) Full-transient closure');
legend({'Short buildup follower','Long buildup follower'},'Location','southeast');
xlim([0 12]); ylim([0 1.1*Cslow]); grid on; box on;
text(0.04,0.96,sprintf('Clearance spread = %.2f m',spread), ...
    'Units','normalized','VerticalAlignment','top','FontSize',9, ...
    'BackgroundColor','w','EdgeColor',[0.85 0.85 0.85]);

ax3=subplot(1,3,3); set(ax3,'Position',[0.705 0.23 0.245 0.63]); hold on;
plot([P.tauLo P.tauHi],[2 2],'-','LineWidth',9,'Color',colors.gray);
plot([loExample hiExample],[1 1],'-','LineWidth',9,'Color',colors.teal);
plot(truthExample,1,'o','MarkerSize',7,'MarkerFaceColor',colors.red, ...
     'MarkerEdgeColor',colors.red);
set(gca,'YTick',[1 2],'YTickLabel',{'D3+','D3'});
xlim([0.15 0.86]); ylim([0.55 2.5]); grid on; box on;
xlabel('Ramp buildup duration, \tau (s)');
ylabel('Measurement contract');
title('(c) Compatible support');
text(0.5,2.32,sprintf('\\eta=%.2f m/s^2; clock +/- %.2f s', ...
    etaExample,clockExample),'HorizontalAlignment','center','FontSize',8);
addSyntheticFooter(fig);
exportFig(fig,out,'Fig1_Physical_ambiguity_and_information');

%% 2. Observation-quality and phase-feasibility map, with correct stopping
etaGrid = [0 0.02 0.05 0.10 0.20 0.35];
clockGrid = [0 0.02 0.04 0.06 0.079 0.10 0.262];
Nedge = 1200;
gain = NaN(numel(clockGrid),numel(etaGrid));
net = NaN(size(gain)); % optional-stop net value; NaN where phase is infeasible
retention = NaN(size(gain));
valid = false(size(gain));
edgeRecords = cell(numel(gain),6);
r=0;
% Common random numbers across assumptions reduce Monte Carlo noise.
rng(20261116, 'twister');
Utruth=rand(Nedge,1); Uno=2*rand(Nedge,1)-1; Uclk=2*rand(Nedge,1)-1;
for iw=1:numel(clockGrid)
    for ie=1:numel(etaGrid)
        eta=etaGrid(ie); w=clockGrid(iw);
        valid(iw,ie) = phaseEligible(P.s,P.p,P.tauLo,P.te,w);
        if valid(iw,ie)
            g=zeros(Nedge,1); good=false(Nedge,1);
            for n=1:Nedge
                truth=P.tauLo+(P.tauHi-P.tauLo)*Utruth(n);
                observed=capFeature(P.v,P.a,P.p,P.s,truth,P.te+w*Uclk(n)) ...
                           + eta*Uno(n);
                [lo,hi,consistent] = posteriorTau(P.a,P.p,P.s, ...
                      P.tauLo,P.tauHi,P.te,eta,w,observed);
                good(n)=consistent && truth>=lo-1e-9 && truth<=hi+1e-9;
                g(n)=max(0,fixedTimeD(P.v,P.a,P.p,P.s,P.tauLo) ...
                         -fixedTimeD(P.v,P.a,P.p,P.s,lo));
            end
            assert(all(good),'Truth retention failed in the generated D3+ data.');
            gain(iw,ie)=mean(g);
            net(iw,ie)=max(0,gain(iw,ie)-P.fee);
            retention(iw,ie)=mean(good);
        end
        r=r+1;
        edgeRecords(r,:)={eta,w,valid(iw,ie),gain(iw,ie), ...
                   net(iw,ie),retention(iw,ie)};
    end
end
Tedge=cell2table(edgeRecords,'VariableNames', ...
    {'sensor_halfwidth_m_s2','clock_halfwidth_s','phase_eligible', ...
     'expected_gross_gap_gain_m','optimal_net_gain_mEquivalent','truth_retained'});
writetable(Tedge,fullfile(out,'Fig2_D3plus_phase_sensor_grid.csv'));

fig=figure('Color','w','Position',[60 80 1650 650]);
ax1=subplot(1,3,1); set(ax1,'Position',[0.055 0.23 0.255 0.63]);
imagesc(etaGrid,1:numel(clockGrid),gain,'AlphaData',double(valid));
set(gca,'YDir','normal','Color',[0.90 0.90 0.90], ...
    'YTick',1:numel(clockGrid),'YTickLabel',cellstr(num2str(clockGrid(:),'%.3f')));
ylabel('Clock half-width (s)'); xlabel('Feature error bound (m/s^2)');
title('(a) Gross information value'); colorbar; box on;

ax2=subplot(1,3,2); set(ax2,'Position',[0.375 0.23 0.255 0.63]);
imagesc(etaGrid,1:numel(clockGrid),net,'AlphaData',double(valid));
set(gca,'YDir','normal','Color',[0.90 0.90 0.90], ...
    'YTick',1:numel(clockGrid),'YTickLabel',cellstr(num2str(clockGrid(:),'%.3f')));
ylabel('Clock half-width (s)'); xlabel('Feature error bound (m/s^2)');
title(sprintf('(b) Net value; fee %.2f m-eq.',P.fee)); colorbar; box on;

ax3=subplot(1,3,3); set(ax3,'Position',[0.705 0.23 0.245 0.63]); hold on;
plot(etaGrid,gain(1,:),'-o','Color',colors.blue,'LineWidth',2.1);
plot(etaGrid,gain(3,:),'-s','Color',colors.teal,'LineWidth',2.1);
plot(etaGrid,gain(5,:),'-d','Color',colors.gray,'LineWidth',2.1);
yline(P.fee,'--','Query fee','Color',colors.red,'LineWidth',1.5);
legend({'0 ms','40 ms','79 ms','Fee'},'Location','northeast');
xlabel('Feature error bound (m/s^2)');
ylabel('Expected gross gain (m)');
title('(c) Stop when value is too low');
grid on; box on;
addSyntheticFooter(fig);
exportFig(fig,out,'Fig2_Phase_feasibility_and_value');

%% 3. A fully calculated three-truck one-query decision reference
caseNames={'Constructed overlapping-response case', ...
           'Synthetic composition of fitted anchors from one truck', ...
           'Narrow-support negative control'};
models=cell(3,1);
models{1}=fleetModel([3.90 3.96 4.02],[.75 .73 .71], ...
                       [.53 .53 .53],[.78 .78 .78]);
% The middle group uses three previously reported configuration-specific
% fitted values from ONE truck in three trailer configurations.
% Composing them into a synthetic fleet is a MODEL ILLUSTRATION ONLY.
models{2}=fleetModel([3.711807 3.836969 4.457561], ...
     [.881604 .615461 .643483], ...
     [.795372 .530793 .530794],[.979560 .683846 .714981]);
models{3}=fleetModel([3.90 3.96 4.02],[.75 .73 .71], ...
                       [.64 .64 .64],[.655 .655 .655]);

feeFleet=P.fee; etaFleet=0.05; clockFleet=0.02;
Nplan=450; Ntest=1100;
allResults=cell(3,1);
policyRows={}; row=0;
for ic=1:3
    model=models{ic};
    fprintf('\nPrecomputing model %d of 3: %s\n',ic,caseNames{ic});
    model=precomputeFleet(model,51,0.010,14.0);
    [baseG,baseOrder,baseGaps]=bestFormation(model, ...
          model.tauLo,model.tauHi,[],0,1);
    model.baseOrder=baseOrder; model.baseG=baseG; model.baseGaps=baseGaps;
    % Separate fitting/evaluation random draws to avoid selecting a query
    % using the same outcomes subsequently presented as validation data.
    planDraws=drawsForFleet(Nplan,70100+ic);
    testDraws=drawsForFleet(Ntest,80100+ic);
    planGain=zeros(1,3); testGain=zeros(1,3);
    testChange=zeros(1,3); idealPlan=zeros(1,3); idealTest=zeros(1,3);
    for q=1:3
        [g1,~,~]=queryGain(model,q,etaFleet,clockFleet,planDraws,false);
        [g2,c2,~]=queryGain(model,q,etaFleet,clockFleet,testDraws,false);
        [gi1,~,~]=queryGain(model,q,etaFleet,clockFleet,planDraws,true);
        [gi2,~,~]=queryGain(model,q,etaFleet,clockFleet,testDraws,true);
        planGain(q)=g1; testGain(q)=g2; testChange(q)=c2;
        idealPlan(q)=gi1; idealTest(q)=gi2;
    end
    [~,selected]=max(planGain);
    [~,selectedIdeal]=max(idealPlan);
    varQ=mostVariableStopDistance(model);
    % Every comparison allows the same one-query budget and fee; stopping
    % is allowed for each executable policy if its EXPECTED gain < fee.
    baseNet=baseG;
    % STOP / TEST actions are selected from PLAN expectations only.
    % Evaluation uses independent draws, including any unfavorable outcomes.
    % Exact expectation over a UNIFORMLY randomly selected truck. Each
    % candidate uses its OWN planned gross gain to decide whether to test.
    % Do not allow an unprofitable chosen truck to borrow value from others.
    randomExpectedGain = mean((planGain > feeFleet).*(testGain-feeFleet));
    randomNet=baseG-randomExpectedGain;
    if planGain(varQ)>feeFleet
        varianceNet=baseG-(testGain(varQ)-feeFleet);
    else
        varianceNet=baseG;
    end
    if planGain(selected)>feeFleet
        selectedNet=baseG-(testGain(selected)-feeFleet);
    else
        selectedNet=baseG;
    end
    if idealPlan(selectedIdeal)>feeFleet
        idealNet=baseG-(idealTest(selectedIdeal)-feeFleet);
    else
        idealNet=baseG;
    end
    allResults{ic}=struct('model',model,'base',baseG,'baseOrder',baseOrder, ...
        'baseGaps',baseGaps,'planGain',planGain,'testGain',testGain, ...
        'testChange',testChange,'selected',selected,'varQ',varQ, ...
        'objectives',[baseNet randomNet varianceNet selectedNet idealNet], ...
        'idealGain',idealTest,'fee',feeFleet);
    labels={'No query / D3','Random (expected; fair stop)','Variance-first D3+', ...
            'Formation-aware D3+','Ideal tau revelation within ramp family (non-deployable)'};
    for k=1:5
        row=row+1;
        policyRows(row,:)={caseNames{ic},labels{k}, ...
             allResults{ic}.objectives(k),baseNet-allResults{ic}.objectives(k), ...
             selected,varQ,feeFleet,etaFleet,clockFleet};
    end
    fprintf('Base order [%s], base clearance %.3f m; test gains %s m\n', ...
        num2str(baseOrder),baseG,mat2str(testGain,3));
end
Tfleet=cell2table(policyRows,'VariableNames', ...
    {'scenario','policy','expected_net_objective_mEquivalent', ...
     'net_gain_vs_no_query_mEquivalent','planned_value_query', ...
     'variance_first_query','one_test_fee_mEquivalent', ...
     'feature_error_halfwidth_m_s2','clock_halfwidth_s'});
writetable(Tfleet,fullfile(out,'Fig3_single_seed_sanity_check.csv'));


% Preserve the predeclared ideal-tau endpoint branch as a supplementary
% record, but do not use it as the main policy-performance illustration.
show=allResults{1}; qDisplay=3;
loCase=show.model.tauLo; hiCase=show.model.tauHi;
hiCase(qDisplay)=loCase(qDisplay);
[branchG,branchOrder,branchGaps]=bestFormation(show.model, ...
    loCase,hiCase,show.baseOrder,show.model.beta,show.base);
branchTable=table(show.base,branchG,show.base-branchG, ...
    {mat2str(show.baseOrder)},{mat2str(branchOrder)}, ...
    {mat2str(show.baseGaps,4)},{mat2str(branchGaps,4)}, ...
    'VariableNames',{'D3_objective_m','IdealTauBranch_objective_m', ...
    'Gross_objective_reduction_m','D3_order', ...
    'IdealTauBranch_order','D3_edge_gaps_m','IdealTauBranch_edge_gaps_m'});
writetable(branchTable,fullfile(out,'Supplementary_predeclared_ideal_tau_branch.csv'));

% Repeat the full planning/evaluation process across independent seeds.
% Confidence intervals below describe Monte Carlo variability under the
% declared synthetic model; they are NOT real-fleet population intervals.
Nrep=30; NplanRep=250; NtestRep=500;
repSeedBase=930000;
R=runRepeatedPolicyValidation(allResults,caseNames,feeFleet,etaFleet, ...
    clockFleet,Nrep,NplanRep,NtestRep,repSeedBase);

writetable(R.repTable,fullfile(out,'Fig3_policy_replications.csv'));
writetable(R.summaryTable,fullfile(out,'Fig3_policy_summary.csv'));
writetable(R.pairedTable,fullfile(out,'Fig3_paired_advantage.csv'));
writetable(R.selectionTable,fullfile(out,'Fig3_query_selection_frequency.csv'));
writetable(R.orderTable,fullfile(out,'Fig3_order_change_summary.csv'));
writetable(R.edgeTable,fullfile(out,'Fig3_per_edge_clearance.csv'));
writetable(R.candidateTable,fullfile(out,'Fig3_candidate_value_audit.csv'));
writetable(R.costRepTable,fullfile(out,'Fig4_cost_frontier_replications.csv'));
writetable(R.costSummaryTable,fullfile(out,'Fig4_cost_frontier_summary.csv'));
writeV4ScientificAudit(R,out,feeFleet,etaFleet,clockFleet);

plotRepeatedPolicyFigure(R,out,colors);
plotRepeatedCostFigure(R,out,colors,feeFleet);



%% 4. Repeated-seed cost frontier and ideal-tau upper benchmark
% The multi-seed cost frontier and nominal-fee comparison are generated by
% plotRepeatedCostFigure. The ideal reference reveals tau exactly within
% the specified one-dimensional ramp family; it is not a general sensor.


%% Final simple numerical checks and provenance
% A policy chosen using finite planning samples may perform worse than
% stopping on independent evaluation draws. Do NOT clip such losses.
assert(all(isfinite(Tfleet.expected_net_objective_mEquivalent)), ...
     'Nonfinite decision objective was generated.');
assert(all(Tedge.truth_retained(~isnan(Tedge.truth_retained))>0.999), ...
     'Physical truth must stay in the assumed support after a valid report.');

fid=fopen(fullfile(out,'DESIGN_NOTES.txt'),'w');
fprintf(fid,'SYNTHETIC MODEL REFERENCE - NOT MEASURED DATA\n\n');
fprintf(fid,'Main theorem family: v=25, a=4, t90=0.9, tau=[0.2,0.8], p=0.9.\n');
fprintf(fid,'Exact D3 residual spread = %.9f m (6.1 m witness).\n',spread);
fprintf(fid,'D3+ sample t_e=0.81 s, phase qualification is explicit.\n');
fprintf(fid,'Measurement error and clock widths are hypothetical settings.\n');
fprintf(fid,'Clock offset 0.262 s can invalidate the common interior-phase contract.\n');
fprintf(fid,'The nominal information fee %.2f m-eq. is not a measured price.\n',P.fee);
fprintf(fid,'Three 3-truck cases: constructed overlapping, fitted-anchor proxies, and narrow control.\n');
fprintf(fid,'The fitted anchors came from one truck in three trailer configurations.\n');
fprintf(fid,'The combined fleet is entirely synthetic, not real simultaneous trucks.\n');
fprintf(fid,'Full-transient edges use ramp displacement at a 0.010-s grid.\n');
fprintf(fid,'Uncertainty extrema in tau are endpoint-valid for p>=0.5 at fixed a and t90.\n');
fprintf(fid,'The grid produces an engineering illustration, not a continuous safety certificate.\n');
fprintf(fid,'The gross query benefit is expected reduction in modeled clearance, not field safety evidence.\n');
fprintf(fid,'The supplementary branch queries truck 3 at its late-onset, shorter-buildup low-tau endpoint.\n');
fprintf(fid,'It is a constructed branch, not the policy-selected or population-representative outcome.\n');
fprintf(fid,'Policy actions (including stopping) are set on separate planning draws.\n');
fprintf(fid,'Independent evaluation can show negative net gains; they are retained.\n');
fprintf(fid,'All deployable policies except forced-variance use optional stopping.\n');
fprintf(fid,'The forced-variance baseline is an ablation with exactly one paid query.\n');
fprintf(fid,'No observed sensor noise law or population frequency is asserted.\n');
fprintf(fid,'No negative/zero cases are filtered out.\n');
fprintf(fid,'V4 multi-seed policy audit: %d independent random replications, %d planning draws and %d test draws per replication.\n',Nrep,NplanRep,NtestRep);
fprintf(fid,'The replications resample outcomes from three fixed synthetic scenarios; they are not independent physical fleets.\n');
fprintf(fid,'Approximate 95%% Student-t intervals quantify Monte Carlo variability only, not real-fleet population uncertainty.\n');
fprintf(fid,'Paired policy differences use common test draws within each replication.\n');
fprintf(fid,'Uniform random chooses a candidate first and stops by THAT candidate plan value.\n');
fprintf(fid,'Expected-random baseline analytically averages across randomized query identities.\n');
fprintf(fid,'Variance-first forced query is an ablation and may have negative net value.\n');
fprintf(fid,'10,000 paired bootstrap resamples yield additional Monte Carlo mean intervals.\n');
fprintf(fid,'Per-edge gap reduction is defined for outcomes retaining the baseline order.\n');
fprintf(fid,'Ideal information reveals exact tau only within the declared ramp family; it is not a deployable full-response sensor.\n');
fprintf(fid,'Scenario 2 combines configuration anchors from one vehicle and is not a simultaneous multi-truck measurement.\n');
fprintf(fid,'Scenario 1 is a deliberately constructed overlap case; scenario 3 is a narrow-support negative control.\n');
fprintf(fid,'The predeclared ideal-tau endpoint branch is stored separately and is not the policy-selected main result.\n');
fprintf(fid,'The cost frontier uses planning-selected queries and independent test gains; negative held-out net gains are retained.\n');
fprintf(fid,'Recommended manuscript use: only as clearly labelled synthetic model illustration.\n');
fclose(fid);

fprintf('\nDone. Reference figures and CSV data saved in:\n%s\n',out);
fprintf('Do NOT present these synthetic figures as empirical validation.\n');
end

% ----------------------------- Local functions --------------------------
function m = fleetModel(a,s,lo,hi)
m = struct('v',25,'p',0.90,'a',a(:)','s',s(:)', ...
   'tauLo',lo(:)','tauHi',hi(:)','delta',0.04,'d0',2.0, ...
   'leadA',6.1,'leadTau',0.38,'leadLag',0.10,'beta',0.012);

% The fitted-anchor example is printed to six decimal places. For its
% second configuration, s - 0.9*tauHi = -4e-7 s solely because of
% decimal rounding. Correct ONLY this sub-microsecond boundary residue;
% reject any materially negative onset instead of relaxing the physics.
roundTol = 1e-6;
onsetAtUpper = m.s - m.p.*m.tauHi;
if any(onsetAtUpper < -roundTol)
    error('Fleet prior has materially negative ramp onset: check s, p and tauHi.');
end
roundingOnly = onsetAtUpper < 0;
if any(roundingOnly)
    originalHi = m.tauHi(roundingOnly);
    m.tauHi(roundingOnly) = (m.s(roundingOnly)-1e-9)./m.p;
    fprintf(['  Corrected tauHi for decimal-rounding onset residue: ', ...
             'max adjustment %.3g s.\n'], max(originalHi-m.tauHi(roundingOnly)));
end
if any(m.tauLo >= m.tauHi)
    error('Fleet prior has empty ramp-duration interval.');
end
if any(m.s - m.p.*m.tauHi < -1e-10)
    error('Fleet prior still has negative ramp onset after rounding correction.');
end
end

function D = fixedTimeD(v,a,p,s,tau)
D = v*s+v.^2./(2*a)+v*(0.5-p).*tau-a.*tau.^2/24;
end

function b = capFeature(~,a,p,s,tau,t)
% The first argument is unused, matching a shared figure-generation call.
lag=s-p.*tau;
b=a.*min(1,max(0,(t-lag)./tau));
end

function x = rampX(t,v,a,tau,lambda,cmd)
% Continuous-time exact position for a delayed, ramp-limited brake profile.
tShift=lambda+cmd;
u=t-tShift;
uStop=v/a+tau/2;
x=v*t;
ramp=(u>0)&(u<=tau);
x(ramp)=v*t(ramp)-a*u(ramp).^3/(6*tau);
plateau=(u>tau)&(u<uStop);
x(plateau)=v*t(plateau)-a/2*(u(plateau)-tau/2).^2-a*tau^2/24;
stopped=(u>=uStop);
stopPos=v*(tShift+uStop)-a/2*(uStop-tau/2)^2-a*tau^2/24;
x(stopped)=stopPos;
end

function yes=phaseEligible(s,p,tauMin,te,clockW)
% Every possible clock-shifted observation lies inside EVERY prior ramp.
maxLag=s-p*tauMin;
yes=(te-clockW>maxLag+1e-12)&&(te+clockW<s-1e-12);
end

function [lo,hi,ok]=posteriorTau(a,p,s,tauMin,tauMax,te,eta,w,y)
% Exact interval inversion for a fixed cap/t90 ramp family in the common
% interior phase. Includes existential bounded noise and clock offset.
lo=tauMin; hi=tauMax; ok=false;
if ~phaseEligible(s,p,tauMin,te,w), return; end
D=s-te;
denomLower=p-(y-eta)/a;
denomUpper=p-(y+eta)/a;
if denomLower<=0, return; end
lo=max(tauMin,(D-w)/denomLower);
if denomUpper>0
    hi=min(tauMax,(D+w)/denomUpper);
end
if lo>hi+1e-9, return; end
lo=min(hi,max(tauMin,lo));
hi=max(lo,min(tauMax,hi));
ok=true;
end

function model=precomputeFleet(model,gridN,dt,tmax)
% A time-grid approximation to the worst full-transient edge closure.
% Endpoint tau choices are physically correct for p>=0.5, fixed cap/t90.
t=0:dt:tmax;
n=3;
model.tauGrid=cell(1,n+1);
X0=cell(1,n+1); XF=cell(1,n+1);
model.tauGrid{1}=model.leadTau;
X0{1}=rampX(t,model.v,model.leadA,model.leadTau,model.leadLag,0)';
XF{1}=X0{1};
for k=1:n
    model.tauGrid{k+1}=linspace(model.tauLo(k),model.tauHi(k),gridN);
    X0{k+1}=zeros(numel(t),gridN);XF{k+1}=zeros(numel(t),gridN);
    for j=1:gridN
       ta=model.tauGrid{k+1}(j);
       lag=model.s(k)-model.p*ta;
       X0{k+1}(:,j)=rampX(t,model.v,model.a(k),ta,lag,0)';
       XF{k+1}(:,j)=rampX(t,model.v,model.a(k),ta,lag,model.delta)';
    end
end
model.C=cell(n+1,n+1);
for i=1:n+1
    for j=2:n+1
        if i==j, continue; end
        Ci=zeros(size(X0{i},2),size(XF{j},2));
        for k=1:size(X0{i},2)
            diff=bsxfun(@minus,XF{j},X0{i}(:,k));
            Ci(k,:)=max(0,max(diff,[],1));
        end
        model.C{i,j}=Ci;
    end
end
model.nGrid=gridN;
end

function [best,chosen,gaps]=bestFormation(model,lower,upper,baseOrder,beta,baseG)
% The worst leading response occurs at high tau (shorter travel);
% the worst following response at low tau (longer travel).
n=3; idxL=ones(1,4);idxH=ones(1,4);
for k=1:n
    grid=model.tauGrid{k+1};
    step=grid(2)-grid(1);
    idxL(k+1)=max(1,min(numel(grid),1+floor((lower(k)-grid(1))/step+1e-10)));
    idxH(k+1)=max(1,min(numel(grid),1+ceil((upper(k)-grid(1))/step-1e-10)));
end
orders=sortrows(perms(1:n));
best=Inf; chosen=[1 2 3];gaps=zeros(1,n);
for r=1:size(orders,1)
    order=orders(r,:);prev=1;gg=zeros(1,n);
    for k=1:n
        j=order(k)+1;
        E=model.C{prev,j};
        gg(k)=model.d0+E(idxH(prev),idxL(j));
        prev=j;
    end
    obj=sum(gg);
    if ~isempty(baseOrder) && beta>0
        obj=obj+beta*baseG*inversionCount(order,baseOrder)/3;
    end
    if obj<best-1e-9
        best=obj;chosen=order;gaps=gg;
    end
end
end

function z=inversionCount(order,baseOrder)
rank=zeros(1,3);
for k=1:3,rank(baseOrder(k))=k;end
q=rank(order);z=0;
for i=1:3
    for j=i+1:3
        z=z+(q(i)>q(j));
    end
end
end

function draw=drawsForFleet(N,seed)
rng(seed,'twister');
draw=struct('truth',rand(N,3),'noise',2*rand(N,3)-1, ...
            'clock',2*rand(N,3)-1);
end

function [gross,change,meanWidth,meanEdgeDelta]=queryGain(model,q,eta,w,draw,ideal)
N=size(draw.truth,1);
val=zeros(N,1);changes=false(N,1); widths=zeros(N,1);
edgeDelta=NaN(N,3);
for k=1:N
    low=model.tauLo;high=model.tauHi;
    trueTau=low(q)+(high(q)-low(q))*draw.truth(k,q);
    if ideal
        low(q)=trueTau;high(q)=trueTau;
    else
        lagmax=model.s(q)-model.p*low(q);
        te=0.5*(lagmax+model.s(q));
        if phaseEligible(model.s(q),model.p,low(q),te,w)
            observed=capFeature(model.v,model.a(q),model.p,model.s(q), ...
                 trueTau,te+w*draw.clock(k,q))+eta*draw.noise(k,q);
            [l,h,ok]=posteriorTau(model.a(q),model.p,model.s(q), ...
                low(q),high(q),te,eta,w,observed);
            if ~ok || trueTau<l-1e-9 || trueTau>h+1e-9
                error('Support-consistency check failed in a fleet draw.');
            end
            low(q)=l;high(q)=h;
        end
    end
    widths(k)=high(q)-low(q);
    [cost,order,gaps]=bestFormation(model,low,high, ...
        model.baseOrder,model.beta,model.baseG);
    val(k)=cost;
    changes(k)=~isequal(order,model.baseOrder);
    if ~changes(k)
        % Comparable physical edges only when the terminal order is intact.
        edgeDelta(k,:)=model.baseGaps-gaps;
    end
end
gross=model.baseG-mean(val);change=mean(changes);meanWidth=mean(widths);
validEdges=~any(isnan(edgeDelta),2);
if any(validEdges)
    meanEdgeDelta=mean(edgeDelta(validEdges,:),1);
else
    meanEdgeDelta=NaN(1,3);
end
if gross<-1e-7
    error('Information unexpectedly increases the expected robust objective.');
end
gross=max(0,gross);
end


function R=runRepeatedPolicyValidation(allResults,caseNames,fee,eta,w, ...
    nRep,nPlan,nTest,seedBase)
% The 30 independent replications resample outcomes on THREE FROZEN SYNTHETIC
% fleet designs. These are not observations from 30 independently sampled
% vehicles/fleets, and the CIs refer only to Monte Carlo uncertainty.
%
% Random methods use ONE query max and the same hypothetical fee:
%  2: Uniformly sampled truck, query-specific stopping.
%  3: Exact expectation over the uniform random truck, query-specific stop.
% The expected-random baseline is Rao-Blackwellized over the action draw;
% it is the main fair random comparison (sampling noise is not rewarded).
% Variance forced is a deliberate ablation; unlike other executable
% policies it must pay to query even when its planned benefit is < fee.

nScenario=numel(allResults);
policyNames={'No query','Random sampled (per-query stop)', ...
    'Random expected (per-query stop)', 'Variance-first (stop)', ...
    'Variance-first (forced)', 'Decision-aware D3+', ...
    'Ideal tau revelation (non-deployable)'};
nPolicy=numel(policyNames);

netGain=NaN(nRep,nScenario,nPolicy);
selectedTruck=NaN(nRep,nScenario,nPolicy);
executionRate=NaN(nRep,nScenario,nPolicy);
routeChange=NaN(nRep,nScenario,nPolicy);
planValue=NaN(nRep,nScenario,nPolicy);
testGross=NaN(nRep,nScenario,nPolicy);
edgeDraw=NaN(nRep,nScenario,3);
repRows={}; rr=0;
candidateRows={}; ar=0;

feeGrid=linspace(0,3.5,71);
nFee=numel(feeGrid);
costValues=NaN(nRep,nScenario,nFee,2);
costRows={}; cr=0;

for ic=1:nScenario
    model=allResults{ic}.model;
    priorWidth=model.tauHi-model.tauLo;
    qVariance=mostVariableStopDistance(model);
    baseOrder=model.baseOrder;
    for ir=1:nRep
        seedPlan=seedBase+ic*100000+ir*100;
        seedTest=seedPlan+1;
        planDraws=drawsForFleet(nPlan,seedPlan);
        testDraws=drawsForFleet(nTest,seedTest);

        pg=zeros(1,3); tg=zeros(1,3);
        tc=zeros(1,3); tw=zeros(1,3); eg=NaN(3,3);
        ig=zeros(1,3); itg=zeros(1,3);
        itc=zeros(1,3); itw=zeros(1,3);
        for q=1:3
            [pg(q),~,~]=queryGain(model,q,eta,w,planDraws,false);
            [tg(q),tc(q),tw(q),eg(q,:)]=queryGain( ...
                model,q,eta,w,testDraws,false);
            [ig(q),~,~]=queryGain(model,q,eta,w,planDraws,true);
            [itg(q),itc(q),itw(q)]=queryGain( ...
                model,q,eta,w,testDraws,true);
        end
        [bestPlan,qBest]=max(pg);
        [bestIdealPlan,qIdeal]=max(ig);
        assert(all(ig>=pg-1e-7), ...
            'The ideal-tau benchmark cannot give less information than D3+.');

        for q=1:3
            ar=ar+1;
            candidateRows(ar,:)={caseNames{ic},ir,q,pg(q),tg(q), ...
                ig(q),itg(q),tc(q),priorWidth(q)-tw(q), ...
                eg(q,1),eg(q,2),eg(q,3),q==qVariance, ...
                seedPlan,seedTest};
        end

        % Random baseline: candidate FIRST, then its OWN stopping criterion.
        rng(seedPlan+37,'twister');
        qRandomDraw=randi(3);
        randomPlan=pg(qRandomDraw);
        if randomPlan>fee
            qRandom=qRandomDraw;
            randomNet=tg(qRandom)-fee;
            randomRoute=tc(qRandom);
            randomGross=tg(qRandom);
        else
            qRandom=0;randomNet=0;randomRoute=NaN;randomGross=NaN;
        end

        % Stronger fair random baseline: analytic expectation over the
        % uniformly selected query, each with its own query-specific stop.
        rndExecute=pg>fee;
        randomExpectedNet=mean(double(rndExecute).*(tg-fee));
        randomExecutionProbability=mean(double(rndExecute));
        if any(rndExecute)
            randomExpectedRoute=mean(tc(rndExecute));
        else
            randomExpectedRoute=NaN;
        end

        % Highest-stopping-distance-variance query, with rational stopping.
        if pg(qVariance)>fee
            qVarianceStop=qVariance;
            varianceNet=tg(qVariance)-fee;
            varianceRoute=tc(qVariance);
            varianceGross=tg(qVariance);
        else
            qVarianceStop=0;varianceNet=0;
            varianceRoute=NaN;varianceGross=NaN;
        end
        % Forced ablation: same query and fee, but NEVER stop.
        forcedNet=tg(qVariance)-fee;

        % Decision-focused query choice with query-specific stopping.
        if bestPlan>fee
            qValue=qBest;valueNet=tg(qValue)-fee;
            valueRoute=tc(qValue);valueGross=tg(qValue);
            edgeDraw(ir,ic,:)=eg(qValue,:);
        else
            qValue=0;valueNet=0;valueRoute=NaN;valueGross=NaN;
        end

        % Ideal revelation of tau within this one-dimensional ramp family.
        if bestIdealPlan>fee
            qIdealExec=qIdeal;
            idealNet=itg(qIdealExec)-fee;
            idealRoute=itc(qIdealExec);
            idealGross=itg(qIdealExec);
        else
            qIdealExec=0;idealNet=0;idealRoute=NaN;idealGross=NaN;
        end

        gains=[0,randomNet,randomExpectedNet,varianceNet,forcedNet, ...
            valueNet,idealNet];
        qs=[0,qRandom,NaN,qVarianceStop,qVariance,qValue,qIdealExec];
        execRate=[0,double(qRandom>0),randomExecutionProbability, ...
            double(qVarianceStop>0),1,double(qValue>0), ...
            double(qIdealExec>0)];
        routes=[NaN,randomRoute,randomExpectedRoute,varianceRoute, ...
            tc(qVariance),valueRoute,idealRoute];
        plans=[0,randomPlan,mean(pg),pg(qVariance),pg(qVariance), ...
            bestPlan,bestIdealPlan];
        grosses=[NaN,randomGross,NaN,varianceGross,tg(qVariance), ...
            valueGross,idealGross];
        for ip=1:nPolicy
            netGain(ir,ic,ip)=gains(ip);
            selectedTruck(ir,ic,ip)=qs(ip);
            executionRate(ir,ic,ip)=execRate(ip);
            routeChange(ir,ic,ip)=routes(ip);
            planValue(ir,ic,ip)=plans(ip);
            testGross(ir,ic,ip)=grosses(ip);

            if isnan(qs(ip))
                action='RANDOMIZED_EXPECTATION';
            elseif qs(ip)==0
                action='STOP';
            else
                action=sprintf('Truck %d',qs(ip));
            end
            rr=rr+1;
            repRows(rr,:)={caseNames{ic},ir,policyNames{ip}, ...
                qs(ip),action,execRate(ip),plans(ip),grosses(ip), ...
                gains(ip),routes(ip),bestPlan,bestIdealPlan, ...
                qVariance,pg(1),pg(2),pg(3),tg(1),tg(2),tg(3), ...
                eta,w,fee,seedPlan,seedTest};
        end

        % Reoptimise the one-test decision to stop at each hypothetical fee.
        % Test draws are frozen independently of planning draws. An expected
        % net result <0 is NOT censored after a query was planned.
        for jf=1:nFee
            cf=feeGrid(jf);
            if bestPlan>cf, vD3=tg(qBest)-cf; else, vD3=0; end
            if bestIdealPlan>cf, vIdeal=itg(qIdeal)-cf; else, vIdeal=0; end
            costValues(ir,ic,jf,1)=vD3;
            costValues(ir,ic,jf,2)=vIdeal;
            for ip=1:2
                cr=cr+1;
                if ip==1
                    name='D3+ decision-aware'; val=vD3;
                else
                    name='Ideal tau revelation (non-deployable)'; val=vIdeal;
                end
                costRows(cr,:)={caseNames{ic},ir,name,cf,val, ...
                    bestPlan,bestIdealPlan,qBest,qIdeal,seedPlan,seedTest};
            end
        end
    end
end

R.repTable=cell2table(repRows,'VariableNames',{ ...
    'scenario','replication','policy','selected_truck_index', ...
    'selected_action','query_execution_probability', ...
    'planning_expected_gross_m','heldout_gross_gain_m', ...
    'heldout_net_gain_m','heldout_order_change_probability', ...
    'best_planning_gross_m','best_ideal_planning_gross_m', ...
    'variance_first_truck','truck1_planning_gross_m', ...
    'truck2_planning_gross_m','truck3_planning_gross_m', ...
    'truck1_heldout_gross_m','truck2_heldout_gross_m', ...
    'truck3_heldout_gross_m','feature_error_halfwidth_m_s2', ...
    'clock_halfwidth_s','query_fee_mEquivalent','planning_seed','test_seed'});
R.candidateTable=cell2table(candidateRows,'VariableNames',{ ...
    'scenario','replication','candidate_truck','planning_gross_m', ...
    'heldout_gross_m','ideal_tau_planning_gross_m', ...
    'ideal_tau_heldout_gross_m','order_change_probability', ...
    'posterior_tau_width_reduction_s','edge1_gain_m', ...
    'edge2_gain_m','edge3_gain_m','variance_first_candidate', ...
    'planning_seed','test_seed'});

nBoot=10000;
summaryRows={};sr=0;
for ic=1:nScenario
    for ip=1:nPolicy
        x=squeeze(netGain(:,ic,ip));
        [mu,med,sd,lo,hi,pos,zer,neg]=repStatsV4(x);
        [blo,bhi]=pairedBootstrapMeanCI(x,nBoot,500000+100*ic+ip);
        sr=sr+1;
        summaryRows(sr,:)={caseNames{ic},policyNames{ip},nRep, ...
            mu,med,sd,lo,hi,blo,bhi,pos,zer,neg, ...
            mean(executionRate(:,ic,ip))};
    end
end
R.summaryTable=cell2table(summaryRows,'VariableNames',{ ...
    'scenario','policy','n_replications','mean_heldout_net_gain_m', ...
    'median_heldout_net_gain_m','sd_between_replications_m', ...
    't95_lower_m','t95_upper_m','bootstrap95_lower_m', ...
    'bootstrap95_upper_m','positive_replication_fraction', ...
    'zero_replication_fraction','negative_replication_fraction', ...
    'query_execution_fraction'});

pairedRows={};pr=0;
pairedRefs={2,'Decision-aware minus random sampled'; ...
            3,'Decision-aware minus random expected'; ...
            4,'Decision-aware minus variance-first stop'; ...
            5,'Decision-aware minus variance forced'};
for ic=1:nScenario
    val=squeeze(netGain(:,ic,6));
    for k=1:size(pairedRefs,1)
        refIndex=pairedRefs{k,1};
        diff=val-squeeze(netGain(:,ic,refIndex));
        [mu,med,sd,lo,hi,pos,zer,neg]=repStatsV4(diff);
        [blo,bhi]=pairedBootstrapMeanCI(diff,nBoot,600000+100*ic+k);
        pr=pr+1;
        pairedRows(pr,:)={caseNames{ic},pairedRefs{k,2},nRep, ...
            mu,med,sd,lo,hi,blo,bhi,pos,zer,neg};
    end
end
R.pairedTable=cell2table(pairedRows,'VariableNames',{ ...
    'scenario','paired_contrast','n_replications','mean_difference_m', ...
    'median_difference_m','sd_difference_m','t95_lower_m', ...
    't95_upper_m','bootstrap95_lower_m','bootstrap95_upper_m', ...
    'positive_fraction','zero_fraction','negative_fraction'});

selectionRows={};qr=0;
selectionPolicies=[2,4,5,6,7];
actions={'Truck 1','Truck 2','Truck 3','STOP'};
acodes=[1,2,3,0];
for ic=1:nScenario
    for ii=1:numel(selectionPolicies)
        ip=selectionPolicies(ii);
        actual=squeeze(selectedTruck(:,ic,ip));
        for ka=1:4
            count=sum(actual==acodes(ka));
            qr=qr+1;
            selectionRows(qr,:)={caseNames{ic},policyNames{ip}, ...
                actions{ka},count,nRep,count/nRep};
        end
    end
end
R.selectionTable=cell2table(selectionRows,'VariableNames',{ ...
    'scenario','policy','decision_action','count','n_replications', ...
    'frequency'});

orderRows={};orow=0;
for ic=1:nScenario
    ex=squeeze(executionRate(:,ic,6));
    x=squeeze(routeChange(ex==1,ic,6));
    if isempty(x)
        mu=NaN;sd=NaN;lo=NaN;hi=NaN;nUse=0;
    else
        [mu,~,sd,lo,hi]=repStatsV4(x);
        nUse=numel(x);
    end
    orow=orow+1;
    orderRows(orow,:)={caseNames{ic},nUse,mu,sd,lo,hi};
end
R.orderTable=cell2table(orderRows,'VariableNames',{ ...
    'scenario','replications_with_query','mean_order_change_probability', ...
    'sd_between_replications','t95_lower','t95_upper'});

edgeRows={};er=0;
for ic=1:nScenario
    bOrder=allResults{ic}.baseOrder;
    for k=1:3
        vals=squeeze(edgeDraw(:,ic,k));
        good=vals(isfinite(vals));
        if isempty(good)
            mu=NaN;med=NaN;sd=NaN;lo=NaN;hi=NaN;nV=0;
        else
            [mu,med,sd,lo,hi]=repStatsV4(good);
            nV=numel(good);
        end
        if k==1, prev=0; else, prev=bOrder(k-1); end
        er=er+1;
        edgeRows(er,:)={caseNames{ic},k,prev,bOrder(k),nV, ...
            mu,med,sd,lo,hi, ...
            R.orderTable.mean_order_change_probability(ic)};
    end
end
R.edgeTable=cell2table(edgeRows,'VariableNames',{ ...
    'scenario','edge_position','leader_truck','follower_truck', ...
    'replications_with_comparable_order','mean_robust_gap_reduction_m', ...
    'median_robust_gap_reduction_m','sd_between_replications_m', ...
    't95_lower_m','t95_upper_m','mean_reordering_probability'});

R.costFeeGrid=feeGrid;
R.costValues=costValues;
costSummary={};cs=0;
for ic=1:nScenario
    for ip=1:2
        for jf=1:nFee
            x=squeeze(costValues(:,ic,jf,ip));
            [mu,med,sd,lo,hi,pos,zer,neg]=repStatsV4(x);
            if ip==1, name='D3+ decision-aware';
            else, name='Ideal tau revelation (non-deployable)'; end
            cs=cs+1;
            costSummary(cs,:)={caseNames{ic},name,feeGrid(jf),nRep, ...
                mu,med,sd,lo,hi,pos,zer,neg};
        end
    end
end
R.costSummaryTable=cell2table(costSummary,'VariableNames',{ ...
    'scenario','policy','fee_mEquivalent','n_replications', ...
    'mean_heldout_net_gain_m','median_heldout_net_gain_m', ...
    'sd_between_replications_m','t95_lower_m','t95_upper_m', ...
    'positive_replication_fraction','zero_replication_fraction', ...
    'negative_replication_fraction'});
R.costRepTable=cell2table(costRows,'VariableNames',{ ...
    'scenario','replication','policy','fee_mEquivalent', ...
    'heldout_net_gain_m','best_planning_gross_m', ...
    'best_ideal_tau_planning_gross_m','selected_D3plus_truck', ...
    'selected_ideal_tau_truck','planning_seed','test_seed'});
R.netGain=netGain;
R.selectedTruck=selectedTruck;
R.executionRate=executionRate;
R.routeChange=routeChange;
R.policyNames=policyNames;
R.caseNames=caseNames;
R.nRep=nRep;
R.nPlan=nPlan;
R.nTest=nTest;
R.nBoot=nBoot;

% Internal data coherence checks - do not claim a physical certificate.
assert(height(R.repTable)==nRep*nScenario*nPolicy);
assert(height(R.candidateTable)==nRep*nScenario*3);
assert(height(R.pairedTable)==nScenario*size(pairedRefs,1));
assert(all(isfinite(netGain(:))));
end

function [mu,med,sd,lo,hi,pos,zer,neg]=repStatsV4(x)
x=x(:); x=x(isfinite(x));
if isempty(x)
    mu=NaN;med=NaN;sd=NaN;lo=NaN;hi=NaN;
    pos=NaN;zer=NaN;neg=NaN;return;
end
mu=mean(x);med=median(x);
if numel(x)>1
    sd=std(x,0);
    half=tCritical95(numel(x)-1)*sd/sqrt(numel(x));
else
    sd=NaN;half=NaN;
end
lo=mu-half;hi=mu+half;
pos=mean(x>1e-10);zer=mean(abs(x)<=1e-10);neg=mean(x<-1e-10);
end

function [lo,hi]=pairedBootstrapMeanCI(x,nBoot,seed)
% Paired-REPLICATION bootstrap: x is already the per-replication paired
% contrast (or a scalar policy performance vector), never raw outcomes.
x=x(:); n=numel(x);
if n==0, lo=NaN;hi=NaN;return;end
if all(x==x(1)),lo=x(1);hi=x(1);return;end
oldState=rng;
rng(seed,'twister');
idx=randi(n,n,nBoot);
bootMean=sort(mean(x(idx),1));
rng(oldState);
lo=interpQuantileSorted(bootMean,0.025);
hi=interpQuantileSorted(bootMean,0.975);
end

function q=interpQuantileSorted(x,p)
n=numel(x);
if n==1,q=x(1);return;end
v=1+(n-1)*p;
i=floor(v);j=ceil(v);
q=x(i)+(v-i)*(x(j)-x(i));
end

function t=tCritical95(df)
% Base MATLAB t critical-value lookup, exact for df 1-30.
tab=[12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 ...
    2.201 2.179 2.160 2.145 2.131 2.120 2.110 2.101 2.093 2.086 ...
    2.080 2.074 2.069 2.064 2.060 2.056 2.052 2.048 2.045 2.042];
if df<1, t=NaN;
elseif df<=numel(tab), t=tab(df);
else, t=1.96;
end
end

function plotRepeatedPolicyFigure(R,out,colors)
scenes=R.caseNames;
shortNames={'Overlap','Fitted proxies','Narrow'};
policies={'Random expected (per-query stop)', ...
    'Variance-first (stop)','Variance-first (forced)', ...
    'Decision-aware D3+','Ideal tau revelation (non-deployable)'};
shortPolicies={'Random expected','Variance-stop','Variance-forced', ...
    'Decision-aware','Ideal tau'};
palette=[colors.gray; 0.69 0.55 0.34; 0.74 0.40 0.34; ...
    colors.blue; colors.teal];

nS=numel(scenes);nP=numel(policies);
means=NaN(nS,nP);half=NaN(nS,nP);
for is=1:nS
    for ip=1:nP
        row=strcmp(R.summaryTable.scenario,scenes{is}) & ...
            strcmp(R.summaryTable.policy,policies{ip});
        means(is,ip)=R.summaryTable.mean_heldout_net_gain_m(row);
        half(is,ip)=R.summaryTable.t95_upper_m(row)-means(is,ip);
    end
end
fig=figure('Color','w','Position',[40 40 1600 980]);
ax1=axes('Position',[0.07 0.57 0.38 0.32]);hold(ax1,'on');
groupW=0.80;bw=groupW/nP;
barHandles=gobjects(1,nP);
for ip=1:nP
    x=(1:nS)-groupW/2+(ip-0.5)*bw;
    barHandles(ip)=bar(ax1,x(:),means(:,ip),bw*0.85, ...
        'FaceColor',palette(ip,:),'EdgeColor','none');
    errorbar(ax1,x(:),means(:,ip),half(:,ip),'k.', ...
        'LineWidth',1.0,'CapSize',6);
end
plot(ax1,[.5 nS+.5],[0 0],'k:','LineWidth',0.9);
set(ax1,'XTick',1:nS,'XTickLabel',shortNames,'FontSize',9);
ylabel(ax1,'Net clearance gain (m-equivalent)');
title(ax1,'(a) Fair policy comparison, 30 seeds','FontSize',11);
lg1=legend(ax1,barHandles,shortPolicies,'Location','northoutside', ...
    'Orientation','horizontal','NumColumns',5,'FontSize',8);
set(lg1,'Units','normalized','Position',[0.085 0.937 0.83 0.038]);
grid(ax1,'on');box(ax1,'on');

ax2=axes('Position',[0.57 0.57 0.36 0.32]);hold(ax2,'on');
contrasts={'Decision-aware minus random expected', ...
           'Decision-aware minus variance-first stop', ...
           'Decision-aware minus variance forced'};
contrastShort={'vs random','vs variance-stop','vs variance-forced'};
cc=[colors.blue; colors.teal;colors.red];
hDiff=gobjects(1,3);
for k=1:3
    xvals=NaN(nS,1);hvals=NaN(nS,1);
    for is=1:nS
        row=strcmp(R.pairedTable.scenario,scenes{is}) & ...
            strcmp(R.pairedTable.paired_contrast,contrasts{k});
        xvals(is)=R.pairedTable.mean_difference_m(row);
        hvals(is)=R.pairedTable.t95_upper_m(row)-xvals(is);
    end
    xx=(1:nS)+(k-2)*0.17;
    hDiff(k)=errorbar(ax2,xx(:),xvals,hvals,'o','Color',cc(k,:), ...
        'LineWidth',1.5,'MarkerFaceColor',cc(k,:),'CapSize',6);
end
plot(ax2,[.5 nS+.5],[0 0],'k:','LineWidth',0.9);
set(ax2,'XTick',1:nS,'XTickLabel',shortNames,'FontSize',9);
ylabel(ax2,'Paired gain advantage (m-equivalent)');
title(ax2,'(b) Matched baseline differences','FontSize',11);
legend(ax2,hDiff,contrastShort,'Location','northeast','FontSize',8);
grid(ax2,'on');box(ax2,'on');

ax3=axes('Position',[0.07 0.14 0.38 0.31]);hold(ax3,'on');
frequencies=zeros(nS,4);
acts={'Truck 1','Truck 2','Truck 3','STOP'};
for is=1:nS
    for ia=1:4
        row=strcmp(R.selectionTable.scenario,scenes{is}) & ...
            strcmp(R.selectionTable.policy,'Decision-aware D3+') & ...
            strcmp(R.selectionTable.decision_action,acts{ia});
        frequencies(is,ia)=R.selectionTable.frequency(row);
    end
end
bars=bar(ax3,1:nS,frequencies,'stacked');
for ia=1:4
    switch ia
        case 1, c=colors.blue;
        case 2, c=colors.teal;
        case 3, c=[0.76 0.57 0.35];
        otherwise, c=[0.55 0.50 0.61];
    end
    set(bars(ia),'FaceColor',c,'EdgeColor','none');
end
set(ax3,'XTick',1:nS,'XTickLabel',shortNames,'FontSize',9);
ylim(ax3,[0 1]);ylabel(ax3,'Selection frequency');
title(ax3,'(c) Decision-aware action selection','FontSize',11);
legend(ax3,bars,acts,'Location','northeast','FontSize',8);
grid(ax3,'on');box(ax3,'on');

ax4=axes('Position',[0.57 0.14 0.36 0.31]);hold(ax4,'on');
sel=strcmp(R.edgeTable.scenario,scenes{1});
E=R.edgeTable(sel,:);
mu=E.mean_robust_gap_reduction_m;
half=E.t95_upper_m-mu;
bar(ax4,1:3,mu,0.55,'FaceColor',colors.teal,'EdgeColor','none');
errorbar(ax4,1:3,mu,half,'k.','LineWidth',1.2,'CapSize',6);
lab=cell(1,3);
for k=1:3
    lab{k}=sprintf('%d to %d',E.leader_truck(k),E.follower_truck(k));
end
set(ax4,'XTick',1:3,'XTickLabel',lab,'FontSize',9);
ylabel(ax4,'Mean conditional edge-gap reduction (m)');
title(ax4,'(d) Gap savings without reordering','FontSize',11);
orderRate=R.orderTable.mean_order_change_probability(1);
if isfinite(orderRate)
    text(ax4,0.04,0.92,sprintf('Order-change probability: %.1f%%', ...
        100*orderRate),'Units','normalized','FontSize',9, ...
        'BackgroundColor','w');
end
grid(ax4,'on');box(ax4,'on');
addSyntheticFooter(fig);
exportFig(fig,out,'Fig3_Fair_policy_comparison_V4');
end

function plotRepeatedCostFigure(R,out,colors,nominalFee)
scenes=R.caseNames;
shortNames={'Overlap','Fitted proxies','Narrow'};
fees=R.costFeeGrid;
fig=figure('Color','w','Position',[50 50 1600 660]);

ax1=axes('Position',[0.08 0.22 0.38 0.65]);hold(ax1,'on');
cc=[colors.blue; colors.teal; colors.gray];
lineHandles=gobjects(1,numel(scenes));
for is=1:numel(scenes)
    a=squeeze(R.costValues(:,is,:,1));
    mu=mean(a,1);
    lo=mu-tCritical95(R.nRep-1)*std(a,0,1)/sqrt(R.nRep);
    hi=mu+tCritical95(R.nRep-1)*std(a,0,1)/sqrt(R.nRep);
    fill(ax1,[fees fliplr(fees)],[lo fliplr(hi)],cc(is,:), ...
        'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
    lineHandles(is)=plot(ax1,fees,mu,'LineWidth',2.1,'Color',cc(is,:));
end
plot(ax1,[min(fees) max(fees)],[0 0],'k:','HandleVisibility','off');
xline(ax1,nominalFee,'--','Nominal fee','Color',[0.55 0.55 0.55], ...
    'HandleVisibility','off');
xlabel(ax1,'Additional measurement fee (m-equivalent)');
ylabel(ax1,'Held-out net gain (m-equivalent)');
title(ax1,'(a) Cost sensitivity and stopping','FontSize',11);
legend(ax1,lineHandles,shortNames,'Location','northeast','FontSize',9);
grid(ax1,'on');box(ax1,'on');

ax2=axes('Position',[0.58 0.22 0.35 0.65]);hold(ax2,'on');
[~,jNom]=min(abs(fees-nominalFee));
palette=[colors.blue; colors.red];
handles=gobjects(1,2);
for ip=1:2
    vals=NaN(numel(scenes),1);half=NaN(numel(scenes),1);
    for is=1:numel(scenes)
        v=squeeze(R.costValues(:,is,jNom,ip));
        vals(is)=mean(v);
        half(is)=tCritical95(R.nRep-1)*std(v,0)/sqrt(R.nRep);
    end
    x=(1:numel(scenes))-0.36+(ip-0.5)*0.36;
    handles(ip)=bar(ax2,x(:),vals,0.30,'FaceColor',palette(ip,:), ...
        'EdgeColor','none');
    errorbar(ax2,x(:),vals,half,'k.','CapSize',6,'LineWidth',1.0);
end
set(ax2,'XTick',1:numel(scenes),'XTickLabel',shortNames,'FontSize',9);
ylabel(ax2,'Held-out net gain (m-equivalent)');
title(ax2,sprintf('(b) Partial vs ideal information; fee %.2f',nominalFee), ...
    'FontSize',11);
lg2=legend(ax2,handles,{'D3+ partial','Ideal tau (not deployable)'}, ...
    'Location','northoutside','Orientation','horizontal','FontSize',9);
set(lg2,'Units','normalized','Position',[0.58 0.91 0.35 0.048]);
grid(ax2,'on');box(ax2,'on');
addSyntheticFooter(fig);
exportFig(fig,out,'Fig4_Cost_and_information_upper_bound_V4');
end

function writeV4ScientificAudit(R,out,fee,eta,w)
fid=fopen(fullfile(out,'V4_SCIENTIFIC_AUDIT.txt'),'w');
if fid<0, error('Cannot open V4 scientific audit output.');end
fprintf(fid,'V4 FAIR SYNTHETIC MODEL VALIDATION - NOT MEASURED DATA\n');
fprintf(fid,'=====================================================\n');
fprintf(fid,'Replications: %d, planning draws: %d, test draws: %d, ', ...
    R.nRep,R.nPlan,R.nTest);
fprintf(fid,'bootstrap resamples: %d\n',R.nBoot);
fprintf(fid,'Fee = %.4f m equivalent; feature error = +/- %.3f m/s^2; ',fee,eta);
fprintf(fid,'clock half-width = %.3f s.\n\n',w);
for ic=1:numel(R.caseNames)
    name=R.caseNames{ic};
    fprintf(fid,'SCENARIO: %s\n',name);
    rows=strcmp(R.summaryTable.scenario,name);
    part=R.summaryTable(rows,:);
    for k=1:height(part)
        fprintf(fid,'  %-41s mean %+9.5f [t95 %.5f, %.5f], ', ...
            part.policy{k},part.mean_heldout_net_gain_m(k), ...
            part.t95_lower_m(k),part.t95_upper_m(k));
        fprintf(fid,'execute probability %.3f\n',part.query_execution_fraction(k));
    end
    ps=strcmp(R.pairedTable.scenario,name);
    pp=R.pairedTable(ps,:);
    for k=1:height(pp)
        fprintf(fid,'  PAIR %-40s %+9.5f [bootstrap %.5f, %.5f]\n', ...
            pp.paired_contrast{k},pp.mean_difference_m(k), ...
            pp.bootstrap95_lower_m(k),pp.bootstrap95_upper_m(k));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'LIMITS: 3 synthetic fixed models, not independently sampled fleets.\n');
fprintf(fid,'The fitted-anchor case combines different configurations of one real truck.\n');
fprintf(fid,'Ideal tau is a within-ramp benchmark, not a deployable sensor.\n');
fprintf(fid,'Repeated Monte Carlo intervals are conditional on the specified model.\n');
fprintf(fid,'Finite-time transient and outer grid calculations are not safety certificates.\n');
fprintf(fid,'Forced-variance is an ablation and can have negative net gains.\n');
fprintf(fid,'Random expected evaluates the fair uniform draw analytically.\n');
fprintf(fid,'No data are tuned to make one method win; all zeros and losses retained.\n');
fclose(fid);
end

function selected=mostVariableStopDistance(model)
varD=zeros(1,3);
for k=1:3
    ts=linspace(model.tauLo(k),model.tauHi(k),201);
    values=fixedTimeD(model.v,model.a(k),model.p,model.s(k),ts);
    varD(k)=var(values,1);
end
[~,selected]=max(varD);
end

function addSyntheticFooter(fig)
annotation(fig,'textbox',[0.08 0.003 0.84 0.030], ...
    'String','SYNTHETIC PHYSICAL-MODEL ILLUSTRATION  |  NOT MEASURED VEHICLE VALIDATION', ...
    'HorizontalAlignment','center','FontSize',9,'Color',[0.44 0.44 0.44], ...
    'LineStyle','none');
end

function exportFig(fig,out,name)
% Export at the on-screen aspect ratio. This avoids portrait-page PDFs
% with a small landscape figure and excessive white margins.
oldUnits = get(fig,'Units');
set(fig,'Units','pixels');
pos = get(fig,'Position');
wIn = pos(3)/120;
hIn = pos(4)/120;
set(fig,'PaperUnits','inches');
set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[0 0 wIn hIn]);
set(fig,'PaperSize',[wIn hIn]);
try
    print(fig,fullfile(out,[name '.pdf']),'-dpdf','-painters');
catch
    print(fig,fullfile(out,[name '.pdf']),'-dpdf');
end
print(fig,fullfile(out,[name '.png']),'-dpng','-r300');
savefig(fig,fullfile(out,[name '.fig']));
set(fig,'Units',oldUnits);
end
