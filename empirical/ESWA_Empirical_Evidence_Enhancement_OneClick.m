function ESWA_Empirical_Evidence_Enhancement_OneClick()
% ESWA_EMPIRICAL_EVIDENCE_ENHANCEMENT_ONECLICK
% Replots the empirical evidence already reported in the current ESWA
% manuscript and Supplement N. No synthetic or fabricated observations
% are generated. No model is refitted by this script.
%
% Source tables in the supplied manuscript:
%   - Main manuscript Table 8 / Supplement N Table N.1:
%       six Freightliner held-out stops, M1 and M3 event-level results.
%   - Supplement N Table N.6:
%       2%, 5%, and 10% pressure-reference threshold summaries.
%   - Supplement N Table N.8:
%       descriptive near-optimal t90 profile widths.
%   - Supplement N Table N.10:
%       Cambridge prefix-extrapolation distance errors.
%
% Purpose:
%   1) Make the speed-fit versus distance-prediction trade-off visible.
%   2) Show whether the conclusion survives the 2/5/10% pressure-reference
%      definitions.
%   3) Visualise the width of near-optimal t90 profiles without calling
%      those widths confidence intervals.
%   4) Display the Cambridge extrapolation stress case without hiding the
%      large Figure 18a errors or scoring the censored Figure 18b trace.
%
% Important:
%   This is a source-data replot, not a new validation experiment.
%   Values are read at full precision from the frozen empirical CSV archive.
%   They were reconciled with Tables N.1, N.6, N.8 and N.10 on 2026-11-19.
%
% MATLAB: executed with R2022b; other releases NOT VERIFIED. Base MATLAB only.
%
% Run:
%   ESWA_Empirical_Evidence_Enhancement_OneClick
%
% Output:
%   ESWA_EMPIRICAL_EVIDENCE_OUTPUT/
%       Figure_Empirical_PairedTradeoff.pdf/.png/.fig
%       Figure_ReferenceThreshold_Robustness.pdf/.png/.fig
%       Figure_NearOptimal_t90_ProfileWidths.pdf/.png/.fig
%       Figure_Cambridge_ExtrapolationStress.pdf/.png/.fig
%       Source_Freightliner_EventLevel.csv
%       Source_Threshold_Sensitivity.csv
%       Source_NearOptimal_ProfileWidths.csv
%       Source_Cambridge_Extrapolation.csv
%       SCIENTIFIC_SCOPE_AUDIT.txt

clc;
close all;

if verLessThan('matlab','9.8')
    error('Vector PDF export requires exportgraphics (R2020a or later).');
end

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
outDir = fullfile(here,'ESWA_EMPIRICAL_EVIDENCE_OUTPUT');
if exist(outDir,'dir') ~= 7
    mkdir(outDir);
end

% Frozen source values: immutable copies from the empirical author archive.
sourceDir = fullfile(here,'frozen_inputs');
H = readtable(fullfile(sourceDir,'pressure_threshold_holdout.csv'),'TextType','string');
S = readtable(fullfile(sourceDir,'pressure_threshold_summary.csv'),'TextType','string');
P = readtable(fullfile(sourceDir,'t90_profile_summary.csv'),'TextType','string');
C = readtable(fullfile(sourceDir,'cambridge_extended_fits.csv'),'TextType','string');
set(groot,'defaultAxesFontName','Arial','defaultTextFontName','Arial', ...
    'defaultAxesFontSize',8.5,'defaultTextFontSize',8.5, ...
    'defaultAxesLineWidth',0.65,'defaultLineLineWidth',1.15);
blue = [0.12 0.36 0.57]; dark = [0.25 0.29 0.33];
palette = [0.78 0.82 0.85; 0.43 0.55 0.63; blue; 0.16 0.23 0.28];

%% ------------------------------------------------------------------------
% 1. Freightliner held-out event-level data: Supplement N Table N.1
% -------------------------------------------------------------------------
event = {'3-1';'3-2';'6-1';'6-2';'10-1';'10-2'};
eventKeys = strcat("run",replace(string(event),"-","_"));

% First-second speed RMSE (m/s)
M1_firstSecond = sourceValues(H,'event',eventKeys,'delayed_step','first_second_rmse_m_s',5,'');
M3_firstSecond = sourceValues(H,'event',eventKeys,'first_order','first_second_rmse_m_s',5,'');

% Signed distance error e_D = predicted - measured (m)
M1_signedDistance = sourceValues(H,'event',eventKeys,'delayed_step','distance_to_3_error_m',5,'');
M3_signedDistance = sourceValues(H,'event',eventKeys,'first_order','distance_to_3_error_m',5,'');

M1_absDistance = abs(M1_signedDistance);
M3_absDistance = abs(M3_signedDistance);

% Positive means M3 improves relative to M1.
earlyRMSEImprovement = M1_firstSecond - M3_firstSecond;
distanceMAEImprovement = M1_absDistance - M3_absDistance;

eventTable = table(event,M1_firstSecond,M3_firstSecond, ...
    earlyRMSEImprovement,M1_signedDistance,M3_signedDistance, ...
    M1_absDistance,M3_absDistance,distanceMAEImprovement, ...
    'VariableNames',{'RunStop','M1_FirstSecondRMSE_mps', ...
    'M3_FirstSecondRMSE_mps','EarlyRMSE_Improvement_M1minusM3_mps', ...
    'M1_SignedDistanceError_m','M3_SignedDistanceError_m', ...
    'M1_AbsoluteDistanceError_m','M3_AbsoluteDistanceError_m', ...
    'AbsoluteDistanceError_Improvement_M1minusM3_m'});
writetable(eventTable,fullfile(outDir,'Source_Freightliner_EventLevel.csv'));

%% Figure 1: paired event-level speed fit and absolute distance error
fig1 = makeFigure(9.0);
tl=tiledlayout(fig1,1,2,'TileSpacing','compact','Padding','compact');
tl.OuterPosition=[0 .12 1 .88];
metrics={M1_firstSecond,M3_firstSecond;M1_absDistance,M3_absDistance};
labels={'First-second speed RMSE (m/s)','Absolute distance error to 3 m/s (m)'};
titles={'(a) Early-speed prediction','(b) Distance prediction'};
for panel=1:2
    ax=nexttile(tl); hold(ax,'on');
    a=metrics{panel,1}; b=metrics{panel,2};
    for k=1:numel(event)
        plot(ax,[a(k) b(k)],[k k],'-','Color',[.65 .69 .72],'HandleVisibility','off');
    end
    h1=plot(ax,a,1:6,'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',4.5);
    h3=plot(ax,b,1:6,'s','Color',dark,'MarkerFaceColor','w','MarkerSize',4.5);
    set(ax,'YTick',1:6,'YTickLabel',event,'YDir','reverse','YLim',[.5 6.5]);
    xlabel(ax,labels{panel}); ylabel(ax,'Run-stop');title(ax,titles{panel},'FontWeight','normal');
    if panel==1, xlim(ax,[0 .65]); else, xlim(ax,[0 10.5]); end
    grid(ax,'on');box(ax,'off');
end
lg=legend([h1 h3],{'M1 delayed step','M3 first-order'},'Orientation','horizontal','Box','off');
lg.Layout.Tile='south';
note(fig1,{'Six paired stops in three held-out runs of one tractor.', ...
    'Paired lines connect models for the same stop; no independent fleet replication.'});
exportFigure(fig1,outDir,'Figure_Empirical_PairedTradeoff');

%% ------------------------------------------------------------------------
% 2. Pressure-reference threshold sensitivity: Supplement N Table N.6
% -------------------------------------------------------------------------
threshold = [2;5;10];

M1_first = sourceValues(S,'pressure_percent',threshold,'delayed_step','first_second_rmse_mean_m_s',[], '');
M3_first = sourceValues(S,'pressure_percent',threshold,'first_order','first_second_rmse_mean_m_s',[], '');

M1_full = sourceValues(S,'pressure_percent',threshold,'delayed_step','speed_rmse_mean_m_s',[], '');
M3_full = sourceValues(S,'pressure_percent',threshold,'first_order','speed_rmse_mean_m_s',[], '');

M1_MAE = sourceValues(S,'pressure_percent',threshold,'delayed_step','distance_mae_m',[], '');
M3_MAE = sourceValues(S,'pressure_percent',threshold,'first_order','distance_mae_m',[], '');

M1_worstUnder = sourceValues(S,'pressure_percent',threshold,'delayed_step','worst_distance_underprediction_m',[], '');
M3_worstUnder = sourceValues(S,'pressure_percent',threshold,'first_order','worst_distance_underprediction_m',[], '');

thresholdTable = table(threshold,M1_first,M3_first,M1_full,M3_full, ...
    M1_MAE,M3_MAE,M1_worstUnder,M3_worstUnder, ...
    M1_first-M3_first,M1_MAE-M3_MAE,M3_worstUnder-M1_worstUnder, ...
    'VariableNames',{'PressureThreshold_percent', ...
    'M1_FirstSecondRMSE_mps','M3_FirstSecondRMSE_mps', ...
    'M1_FullWindowRMSE_mps','M3_FullWindowRMSE_mps', ...
    'M1_DistanceMAE_m','M3_DistanceMAE_m', ...
    'M1_LargestUnderprediction_m','M3_LargestUnderprediction_m', ...
    'EarlyRMSE_Improvement_M1minusM3_mps', ...
    'DistanceMAE_Improvement_M1minusM3_m', ...
    'WorstUnderprediction_Increase_M3minusM1_m'});
writetable(thresholdTable,fullfile(outDir,'Source_Threshold_Sensitivity.csv'));

fig2 = makeFigure(8.4);
tl=tiledlayout(fig2,1,3,'TileSpacing','compact','Padding','compact');
tl.OuterPosition=[0 .15 1 .85];
metrics={M1_first,M3_first;M1_MAE,M3_MAE;M1_worstUnder,M3_worstUnder};
labels={'First-second RMSE (m/s)','Distance MAE (m)','Largest underprediction (m)'};
titles={'(a) Early-speed prediction','(b) Distance prediction','(c) Worst underprediction'};
for panel=1:3
    ax=nexttile(tl);hold(ax,'on');
    h1=plot(ax,threshold,metrics{panel,1},'-o','Color',blue,'MarkerSize',4.5,'MarkerFaceColor',blue);
    h3=plot(ax,threshold,metrics{panel,2},'--s','Color',dark,'MarkerSize',4.5,'MarkerFaceColor','w');
    set(ax,'XTick',threshold,'XLim',[1 11]);
    xlabel(ax,{'Pressure threshold','(%)'});ylabel(ax,labels{panel});
    title(ax,titles{panel},'FontSize',8.5,'FontWeight','normal');grid(ax,'on');box(ax,'off');
end
lg=legend([h1 h3],{'M1 delayed step','M3 first-order'},'Orientation','horizontal','Box','off');
lg.Layout.Tile='south';
note(fig2,{'Each threshold changes the event origin and fitting interval.', ...
    'Scores are descriptive; distance errors are not pooled across origins.'});
exportFigure(fig2,outDir,'Figure_ReferenceThreshold_Robustness');

%% ------------------------------------------------------------------------
% 3. Near-optimal t90 profile widths: Supplement N Table N.8
% -------------------------------------------------------------------------
trace = {'2-1';'2-2';'3-1';'3-2';'5-1';'5-2';'6-1';'6-2'; ...
    '9-1';'9-2';'10-1';'10-2';'Fig.13a';'Fig.13b'; ...
    'Fig.15a';'Fig.15b';'Fig.18a';'Fig.18b'};
traceKeys = replace(strcat("run",replace(string(trace),"-","_")),"runFig.","figure_");

M2_t90 = sourceValues(P,'event',traceKeys,'ramp','free_t90_s',5,'');
M2_width = sourceValues(P,'event',traceKeys,'ramp','compatible_grid_span_s',5,'');

M3_t90 = sourceValues(P,'event',traceKeys,'first_order','free_t90_s',5,'');
M3_width = sourceValues(P,'event',traceKeys,'first_order','compatible_grid_span_s',5,'');

profileTable = table(trace,M2_t90,M2_width,M3_t90,M3_width, ...
    'VariableNames',{'Trace','M2_Ramp_t90_s','M2_NearOptimal_t90Span_s', ...
    'M3_FirstOrder_t90_s','M3_NearOptimal_t90Span_s'});
writetable(profileTable,fullfile(outDir,'Source_NearOptimal_ProfileWidths.csv'));

fig3 = makeFigure(9.2);
ax=axes(fig3,'Position',[.10 .27 .88 .64]);
b=bar(ax,[M2_width M3_width],'grouped','EdgeColor','none');
b(1).FaceColor=blue;b(2).FaceColor=[.60 .66 .71];
set(ax,'XTick',1:numel(trace),'XTickLabel',trace,'FontSize',8.5,'XLim',[.3 18.7]);
xtickangle(ax,45);ylabel(ax,'Near-optimal t_{90} grid span (s)');
title(ax,{'Timing-profile sensitivity at a 0.02 m/s RMSE tolerance'},'FontWeight','normal');
xline(ax,12.5,':','Color',[.4 .4 .4],'HandleVisibility','off');
legend(ax,{'M2 ramp','M3 first-order'},'Location','northwest','Box','off');
grid(ax,'on');box(ax,'off');
note(fig3,{'Freightliner: 5% pressure reference; Cambridge: demand reference.', ...
    'Descriptive grid spans, not confidence intervals or identified parameter bounds.'});
exportFigure(fig3,outDir,'Figure_NearOptimal_t90_ProfileWidths');

%% ------------------------------------------------------------------------
% 4. Cambridge prefix extrapolation stress: Supplement N Table N.10
% -------------------------------------------------------------------------
camTrace = {'Fig.13a';'Fig.13b';'Fig.15a';'Fig.15b';'Fig.18a'};
camKeys = replace(string(camTrace),"Fig.","figure_");
censored = C.fit_window=="prefix" & C.event=="figure_18b";
assert(nnz(censored)==4 && all(isnan(C.distance_to_3_error_m(censored))), ...
    'Censored Figure 18b must remain unscored.');
M0_error = sourceValues(C,'event',camKeys,'constant','distance_to_3_error_m',[],'prefix');
M1_error = sourceValues(C,'event',camKeys,'delayed_step','distance_to_3_error_m',[],'prefix');
M2_error = sourceValues(C,'event',camKeys,'ramp','distance_to_3_error_m',[],'prefix');
M3_error = sourceValues(C,'event',camKeys,'first_order','distance_to_3_error_m',[],'prefix');

camTable = table(camTrace,M0_error,M1_error,M2_error,M3_error, ...
    'VariableNames',{'Trace','M0_Error_m','M1_Error_m','M2_Error_m','M3_Error_m'});
writetable(camTable,fullfile(outDir,'Source_Cambridge_Extrapolation.csv'));

fig4 = makeFigure(9.0);
tl=tiledlayout(fig4,1,2,'TileSpacing','compact','Padding','compact');
tl.OuterPosition=[0 .15 1 .85];
ax=nexttile(tl);b=bar(ax,[M0_error(1:4) M1_error(1:4) M2_error(1:4) M3_error(1:4)],'grouped','EdgeColor','none');
for k=1:4,b(k).FaceColor=palette(k,:);end
set(ax,'XTick',1:4,'XTickLabel',camTrace(1:4),'YLim',[-.6 6.2]);
ylabel(ax,'Signed distance error (m)');title(ax,'(a) Four other scored traces','FontWeight','normal');
yline(ax,0,'-','Color',dark,'HandleVisibility','off');grid(ax,'on');box(ax,'off');
lg=legend(ax,{'M0','M1','M2','M3'},'Orientation','horizontal','Box','off');lg.Layout.Tile='south';
ax=nexttile(tl);vals=[M0_error(5) M1_error(5) M2_error(5) M3_error(5)];
b=bar(ax,vals,'FaceColor','flat','EdgeColor','none');b.CData=palette;
set(ax,'XTick',1:4,'XTickLabel',{'M0','M1','M2','M3'},'YLim',[0 110]);
ylabel(ax,'Signed distance error (m)');title(ax,'(b) Fig. 18a stress case','FontWeight','normal');
for k=1:4,text(ax,k,vals(k)+3,sprintf('%.3f',vals(k)),'HorizontalAlignment','center','FontSize',8);end
grid(ax,'on');box(ax,'off');
note(fig4,{'Positive = predicted minus measured distance; all Fig. 18a errors exceed +93 m.', ...
    'Fig. 18b is censored before 3 m/s and has no scored distance error.'});
exportFigure(fig4,outDir,'Figure_Cambridge_ExtrapolationStress');

%% ------------------------------------------------------------------------
% Audit file
% -------------------------------------------------------------------------
fid = fopen(fullfile(outDir,'SCIENTIFIC_SCOPE_AUDIT.txt'),'w');
fprintf(fid,'ESWA empirical evidence replot audit\n');
fprintf(fid,'====================================\n\n');
fprintf(fid,'No new observations, model fits, Monte Carlo outcomes or synthetic data were generated.\n');
fprintf(fid,'All values are read from frozen empirical CSVs reconciled with tables N.1, N.6, N.8 and N.10.\n');
fprintf(fid,'The accompanying DATA_RECONCILIATION_REPORT.csv documents every original transcription.\n\n');
fprintf(fid,'Freightliner event-level interpretation:\n');
fprintf(fid,'- First-second RMSE improves in %d of 6 held-out stops.\n',sum(earlyRMSEImprovement>1e-12));
fprintf(fid,'- Absolute distance error improves in %d of 6 held-out stops.\n',sum(distanceMAEImprovement>1e-12));
fprintf(fid,'- The 6 events are clustered in 3 held-out runs; do not treat them as 6 independent vehicles.\n\n');
fprintf(fid,'Threshold sensitivity:\n');
fprintf(fid,'- M3 first-second RMSE improvement over M1: %.3f, %.3f, %.3f m/s at 2, 5, 10%%.\n', ...
    M1_first(1)-M3_first(1),M1_first(2)-M3_first(2),M1_first(3)-M3_first(3));
fprintf(fid,'- M3 distance MAE improvement over M1: %.3f, %.3f, %.3f m at 2, 5, 10%%.\n', ...
    M1_MAE(1)-M3_MAE(1),M1_MAE(2)-M3_MAE(2),M1_MAE(3)-M3_MAE(3));
fprintf(fid,'- M3 increase in largest underprediction: %.3f, %.3f, %.3f m at 2, 5, 10%%.\n', ...
    M3_worstUnder(1)-M1_worstUnder(1), ...
    M3_worstUnder(2)-M1_worstUnder(2), ...
    M3_worstUnder(3)-M1_worstUnder(3));
fprintf(fid,'- Near-optimal t90 spans are descriptive grid spans, not confidence intervals.\n');
fprintf(fid,'- Cambridge Fig. 18a is a retrospective extrapolation stress case; Fig. 18b has no scored endpoint.\n');
fprintf(fid,'- Do not present these re-plots as new independent validation.\n');
fclose(fid);

fprintf('\nCompleted source-data replot.\n');
fprintf('Output folder: %s\n',outDir);
fprintf('No new experiments or fabricated outcomes were generated.\n\n');

end

function values = sourceValues(T,keyName,keys,model,fieldName,pressure,window)
values=zeros(numel(keys),1);
for k=1:numel(keys)
    if isnumeric(keys), mask=T.(keyName)==keys(k); else, mask=T.(keyName)==string(keys(k)); end
    mask=mask & T.model==string(model);
    if ~isempty(pressure) && ismember('pressure_percent',T.Properties.VariableNames)
        % Cambridge has a demand clock, not a Freightliner pressure percentage.
        mask=mask & (T.pressure_percent==pressure | T.dataset=="Cambridge");
    end
    if ~isempty(window),mask=mask & T.fit_window==string(window);end
    assert(nnz(mask)==1,'Source row is missing or duplicated.');
    values(k)=T.(fieldName)(mask);
    assert(isfinite(values(k)),'A scored field must be finite.');
end
end

function fig=makeFigure(height)
fig=figure('Color','w','Visible','off','Units','centimeters','Position',[2 2 16.6 height]);
end

function note(fig,lines)
annotation(fig,'textbox',[.04 .005 .94 .09],'String',lines,'EdgeColor','none', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8,'Interpreter','none');
end

function exportFigure(fig,outDir,name)
set(fig,'Color','w','Renderer','painters','PaperUnits','centimeters');
pos=get(fig,'Position');set(fig,'PaperSize',pos(3:4),'PaperPosition',[0 0 pos(3:4)]);
exportgraphics(fig,fullfile(outDir,[name '.pdf']),'ContentType','vector','BackgroundColor','white');
print(fig,fullfile(outDir,[name '.png']),'-dpng','-r300');
savefig(fig,fullfile(outDir,[name '.fig']));
end
