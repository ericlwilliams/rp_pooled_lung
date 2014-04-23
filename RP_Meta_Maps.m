function RP_Meta_Maps
%% %%%%%%%%%%%%%%%%%%%%%%%%
% Plotting combined data likelihood maps, two ways:
% - slice combined llhd matrix by max combined llhd, project meta
% - slice combined llhd matrix at meta points, slice of combined llhd NOT
% projection
%%%%%%%%%%%%%%%%%%%%%%%%

%% random effects meta analysis results [val,95% lo, hi]
re_td50 = [21.88,15.9,27.8];
re_log10a = [-0.0862,-0.69897,0.15533];
re_m = [0.29,0.13,0.44];

%% fixed effects 
fe_td50 = [21.88,15.9,27.8];
fe_log10a = [-0.0506,-0.2924,0.1038];
fe_m = [0.21,0.14,0.28];

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';

% Only loading combined dataset 
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
%fn='C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_med_EUD_fine_meta_comb.mat';
load(fn,'CGcomb');

% LymanN = log10(CGcomb.mLymanN);
% CGcomb.mLymanN = LymanN;


% use mx_llhd_llhd to find contours when slicing at meta points
[mx_llhd,mx_llhd_ind] = max(CGcomb.mLymanGrid.loglikelihood(:));
[comb_td50_idx,comb_m_idx,comb_n_idx] =...
    ind2sub(size(CGcomb.mLymanGrid.loglikelihood),mx_llhd_ind);
% compute the 68% and 95% CIs
% uses chi2inv(0.65,1)...
% low68 = mx_llhd-0.5* (1);
% low95 = mx_llhd-0.5* (1.96^2);
% low99 = mx_llhd-0.5* (3^2);

%double check
low68 = mx_llhd -0.5*(chi2inv(0.68,2));
low95 = mx_llhd -0.5*(chi2inv(0.95,2));
low99 = mx_llhd -0.5*(chi2inv(0.99,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TD50 vs m 
% split at meta analysis n
%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_n = 10.^(-fe_log10a);
%find position of meta n
[~,n_idx] = min(abs(CGcomb.mLymanN-fe_n(1)));   
llhds = CGcomb.mLymanGrid.loglikelihood(:,:,n_idx);


cur_fig=figure(1); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanGrid.TD50,CGcomb.mLymanGrid.m,llhds',[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanGrid.TD50(comb_td50_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_td50(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_td50(1) fe_td50(1)],[fe_m(2) fe_m(3)],'k','LineWidth',1);
plot([fe_td50(2) fe_td50(3)],[fe_m(1) fe_m(1)],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit (projected)','Meta analysis best fit',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([15 28]);
ylim([0.1 0.6]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('TD_{50} (Gy)','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_meta_td50_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_meta_td50_m.png...']);
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TD50 vs a 
% split at meta analysis m
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-fe_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(2); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanN,CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_n(1),fe_td50(1),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_n(1) fe_n(1)],[fe_td50(2) fe_td50(3)],'k','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_td50(1) fe_td50(1)],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit (projected)','Meta analysis best fit',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.7 3]);
ylim([10 33]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_meta_n_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_meta_n_td50.png...']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot m vs a 
% split at meta analysis TD50
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,td50_idx] = min(abs(CGcomb.mLymanGrid.TD50-fe_td50(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(td50_idx,:,:));

cur_fig=figure(3); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_n(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_n(1) fe_n(1)],[fe_m(2) fe_m(3)],'k','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_m(1) fe_m(1)],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit (projected)','Meta analysis best fit',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.7 3.5]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_meta_n_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_meta_n_m.png...']);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%
% slice at best combined model
% project meta analysis model
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% TD50 vs m 
%%%%%%%%%%%%%%%%%%%%%%%%%%

llhds = CGcomb.mLymanGrid.loglikelihood(:,:,comb_n_idx);

cur_fig=figure(10); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanGrid.TD50,CGcomb.mLymanGrid.m,llhds',[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanGrid.TD50(comb_td50_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_td50(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_td50(1) fe_td50(1)],[fe_m(2) fe_m(3)],'k','LineWidth',1);
plot([fe_td50(2) fe_td50(3)],[fe_m(1) fe_m(1)],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit','Meta analysis best fit (projected)',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([15 28]);
ylim([0.1 0.6]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('TD_{50} (Gy)','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_td50_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_td50_m.png...']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot n vs TD50
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,comb_m_idx,:));

cur_fig=figure(20); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanN,CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_n(1),fe_td50(1),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_n(1) fe_n(1)],[fe_td50(2) fe_td50(3)],'k','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_td50(1) fe_td50(1)],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit','Meta analysis best fit (projected)',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.7 3]);
ylim([10 33]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_n_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_n_td50.png...']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot n vs m
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(comb_td50_idx,:,:));

cur_fig=figure(30); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_n(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_n(1) fe_n(1)],[fe_m(2) fe_m(3)],'k','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_m(1) fe_m(1)],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit','Meta analysis best fit (projected)',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.7 3.5]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_n_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_n_m.png...']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot n vs log10(TD50)
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,comb_m_idx,:));

cur_fig=figure(40); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanN,log10(CGcomb.mLymanGrid.TD50),llhds,[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanN(comb_n_idx),log10(CGcomb.mLymanGrid.TD50(comb_td50_idx)),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_n(1),log10(fe_td50(1)),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_n(1) fe_n(1)],[log10(fe_td50(2)) log10(fe_td50(3))],'k','LineWidth',1);
plot([fe_n(2) fe_n(3)],[log10(fe_td50(1)) log10(fe_td50(1))],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit','Meta analysis best fit (projected)',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.7 3]);
ylim([1.05 1.55]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('log(TD50) (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_n_log10_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_n_log10_td50.png...']);
end


%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-fe_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(4); clf reset;
set(cur_fig,'Position',ss_four2three);

contour(CGcomb.mLymanN,log10(CGcomb.mLymanGrid.TD50),llhds,[low99,low95,low68],'LineWidth',2);hold on; 
h_comb=plot(CGcomb.mLymanN(comb_n_idx),log10(CGcomb.mLymanGrid.TD50(comb_td50_idx)),'r+','LineWidth',2,'MarkerSize',12);
h_meta=plot(fe_n(1),log10(fe_td50(1)),'kx','LineWidth',2,'MarkerSize',12);
plot([fe_n(1) fe_n(1)],[log10(fe_td50(2)) log10(fe_td50(3))],'k','LineWidth',1);
plot([fe_n(2) fe_n(3)],[log10(fe_td50(1)) log10(fe_td50(1))],'k','LineWidth',1);

lgnd = legend([h_comb h_meta],'Pooled best fit (projected)','Meta analysis best fit',...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.7 3]);
ylim([1.05 1.55]);

set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('log(TD50) (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_meta_n_log10_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_meta_n_log10_td50.png...']);
end

end