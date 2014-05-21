function RP_Meta_Maps
%% %%%%%%%%%%%%%%%%%%%%%%%%
% Plotting combined data likelihood maps, two ways:
% - slice combined llhd matrix by max combined llhd, project meta
% - slice combined llhd matrix at meta points, slice of combined llhd NOT
% projection
%%%%%%%%%%%%%%%%%%%%%%%%
one_sided_var = true;

cntr_colors = {[1 0 0], [0 1 0], [0 0 1]};
% cm1 = colormap(jet(300)); cm1=cm1(1:256,:); %cm1(end,:) = 0.5;
% cm2 = colormap(jet(10));

%% %%%%%%%%%%%%%%%%%%%%%%%%
% Meta analysis results
% data from RP_Meta_Analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed effects meta analysis results [val,95% CIs]

% two-sided variance
fe_a = [0.8878,0.5052,1.27];
fe_n = [0.7574,0.4394,1.075];
fe_td50 = [21.89,15.96,27.81];
fe_m = [0.2144,0.14,0.2888];

fe_log10a = log10(fe_a);

%% random effects meta analysis results [val,95% CIs]

% two-sided variance
re_a = [0.8169,0.2039,1.43];
re_n = [1.237,0.04703,2.427];
re_td50 = [21.89,15.96,27.81];
re_m = [0.287,0.1311,0.4428];

re_log10a = log10(re_a);


screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';

% Only loading combined dataset 
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
%fn='C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_med_EUD_fine_meta_comb.mat';

if isunix %on mac
    fig_loc = strrep(fig_loc,'Z:/elw/','/Users/elw/Documents/');
    fn = strrep(fn,'Z:/elw/','/Users/elw/Documents/');
end


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
low68 = mx_llhd -0.5*(chi2inv(0.68,3));
low95 = mx_llhd -0.5*(chi2inv(0.95,3));
low99 = mx_llhd -0.5*(chi2inv(0.99,3));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TD50 vs m 
% split at meta analysis n
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,n_idx] = min(abs(CGcomb.mLymanN-fe_n(1)));   
llhds = CGcomb.mLymanGrid.loglikelihood(:,:,n_idx);

cur_fig=figure(1); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(CGcomb.mLymanGrid.TD50,CGcomb.mLymanGrid.m,llhds',[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');

for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanGrid.TD50(comb_td50_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_td50(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_td50(1) fe_td50(1)],[fe_m(2) fe_m(3)],'k.-','LineWidth',1);
plot([fe_td50(2) fe_td50(3)],[fe_m(1) fe_m(1)],'k.-','LineWidth',1);

lgnd = legend([h_comb h_fe_meta],...
    'Pooled best fit (projected)',...
    ['Fixed-effects',10,'meta analysis best fit'],...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([15 32]);
ylim([0.04 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('TD_{50} (Gy)','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_fe_meta_td50_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_fe_meta_td50_m.png...']);
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot n vs TD50
% split at meta analysis m
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-fe_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(2); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(CGcomb.mLymanN,CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_n(1),fe_td50(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_n(1) fe_n(1)],[fe_td50(2) fe_td50(3)],'k.-','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_td50(1) fe_td50(1)],'k.-','LineWidth',1);
lgnd = legend([h_comb h_fe_meta],...
    'Pooled best fit (projected)',...
    ['Fixed-effects',10,'meta analysis best fit'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);


xlim([0.26 3.2]);
ylim([12 34.5]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_fe_meta_n_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_fe_meta_n_td50.png...']);
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

[~,Htmp]=contour(CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_n(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_n(1) fe_n(1)],[fe_m(2) fe_m(3)],'k.-','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_m(1) fe_m(1)],'k.-','LineWidth',1);

lgnd = legend([h_comb h_fe_meta],...
    'Pooled best fit (projected)',...
    ['Fixed-effects',10,'meta analysis best fit'],...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.26 4.4]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_fe_meta_n_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_fe_meta_n_m.png...']);
end

%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-fe_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(4); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(log10(CGcomb.mLymanN),log10(CGcomb.mLymanGrid.TD50),llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(log10(CGcomb.mLymanN(comb_n_idx)),log10(CGcomb.mLymanGrid.TD50(comb_td50_idx)),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(log10(fe_n(1)),log10(fe_td50(1)),'kx','LineWidth',2,'MarkerSize',14);
plot(log10([fe_n(1) fe_n(1)]),[log10(fe_td50(2)) log10(fe_td50(3))],'k.-','LineWidth',1);
plot(log10([fe_n(2) fe_n(3)]),[log10(fe_td50(1)) log10(fe_td50(1))],'k.-','LineWidth',1);

lgnd = legend([h_comb h_fe_meta],...
    'Pooled best fit (projected)',...
    ['Fixed-effects',10,'meta analysis best fit'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([-0.2 0.75]);
ylim([1.1 1.55]);


set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('log_{10}(n)','FontSize',22); ylabel('log(TD50) (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_fe_meta_log10_n_log10_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_fe_meta_log10_n_log10_td50.png...']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a vs TD50
% split at meta analysis m
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-fe_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(5); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour((1./CGcomb.mLymanN),CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(1./(CGcomb.mLymanN(comb_n_idx)),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_a(1),fe_td50(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_a(1) fe_a(1)],[fe_td50(2) fe_td50(3)],'k.-','LineWidth',1);
plot([fe_a(2) fe_a(3)],[fe_td50(1) fe_td50(1)],'k.-','LineWidth',1);
lgnd = legend([h_comb h_fe_meta],...
    'Pooled best fit (projected)',...
    ['Fixed-effects',10,'meta analysis best fit'],...
    'Location','NorthWest');
set(lgnd,'FontSize',18);


xlim([0.1 1.5]);
ylim([12 34.5]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman a','FontSize',22); 
ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_fe_meta_a_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_fe_meta_a_td50.png...']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot m vs a 
% split at meta analysis TD50
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,td50_idx] = min(abs(CGcomb.mLymanGrid.TD50-fe_td50(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(td50_idx,:,:));

cur_fig=figure(6); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(1./CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(1./CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_a(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_a(1) fe_a(1)],[fe_m(2) fe_m(3)],'k.-','LineWidth',1);
plot([fe_a(2) fe_a(3)],[fe_m(1) fe_m(1)],'k.-','LineWidth',1);

lgnd = legend([h_comb h_fe_meta],...
    'Pooled best fit (projected)',...
    ['Fixed-effects',10,'meta analysis best fit'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([0.1 1.5]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman a','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_fe_meta_a_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_fe_meta_a_m.png...']);
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

% TMP
%  low68 = mx_llhd - 0.5* (1);

[Ctmp,Htmp]=contour(CGcomb.mLymanGrid.TD50,CGcomb.mLymanGrid.m,llhds',...
    [low99,low95,low68],'LineWidth',2);hold on;
Cld = get(Htmp,'Children');

% % TMP
% c_td50 = Ctmp(1,:);
% [~,c_68_idx] = min(abs(c_td50-low68));
% [~,c_95_idx] = min(abs(c_td50-low95));
% [~,c_99_idx] = min(abs(c_td50-low99));
% c_td50_low68 = c_td50(c_68_idx+1:end);
% c_td50_low95 = c_td50(c_95_idx+1:c_68_idx-1);
% c_td50_low99 = c_td50(c_99_idx+1:c_95_idx-1);


for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanGrid.TD50(comb_td50_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_td50(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_td50(1) fe_td50(1)],[fe_m(2) fe_m(3)],'k.-.','LineWidth',1);
plot([fe_td50(2) fe_td50(3)],[fe_m(1) fe_m(1)],'k.-.','LineWidth',1);
h_re_meta=plot(re_td50(1),re_m(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_td50(1) re_td50(1)],[re_m(2) re_m(3)],'m.-.','LineWidth',1);
plot([re_td50(2) re_td50(3)],[re_m(1) re_m(1)],'m.-.','LineWidth',1);

lgnd = legend([h_comb h_fe_meta h_re_meta],...
    'Pooled best fit',...
    ['Fixed-effects meta analysis',10,'best fit (projected)'],...
    ['Random-effects meta analysis',10,'best fit (projected)'],...
    'Location','SouthEast');

set(lgnd,'FontSize',18);

xlim([15 32]);
ylim([0.04 0.65]);
 
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

[~,Htmp]=contour(CGcomb.mLymanN,CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_n(1),fe_td50(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_n(1) fe_n(1)],[fe_td50(2) fe_td50(3)],'k.-.','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_td50(1) fe_td50(1)],'k.-.','LineWidth',1);
h_re_meta=plot(re_n(1),re_td50(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_n(1) re_n(1)],[re_td50(2) re_td50(3)],'m.-.','LineWidth',1);
plot([re_n(2) re_n(3)],[re_td50(1) re_td50(1)],'m.-.','LineWidth',1);

lgnd = legend([h_comb h_fe_meta h_re_meta],...
    'Pooled best fit',...
    ['Fixed-effects meta analysis',10,'best fit (projected)'],...
    ['Random-effects meta analysis',10,'best fit (projected)'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([0.26 3.2]);
ylim([12 34.5]);
 
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

[~,Htmp]=contour(CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_n(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_n(1) fe_n(1)],[fe_m(2) fe_m(3)],'k.-.','LineWidth',1);
plot([fe_n(2) fe_n(3)],[fe_m(1) fe_m(1)],'k.-.','LineWidth',1);
h_re_meta=plot(re_n(1),re_m(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_n(1) re_n(1)],[re_m(2) re_m(3)],'m.-.','LineWidth',1);
plot([re_n(2) re_n(3)],[re_m(1) re_m(1)],'m.-.','LineWidth',1);

lgnd = legend([h_comb h_fe_meta h_re_meta],...
    'Pooled best fit',...
    ['Fixed-effects meta analysis',10,'best fit (projected)'],...
    ['Random-effects meta analysis',10,'best fit (projected)'],...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.26 4.4]);
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

[~,Htmp]=contour(log10(CGcomb.mLymanN),log10(CGcomb.mLymanGrid.TD50),llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(log10(CGcomb.mLymanN(comb_n_idx)),log10(CGcomb.mLymanGrid.TD50(comb_td50_idx)),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(-fe_log10a(1),log10(fe_td50(1)),'kx','LineWidth',2,'MarkerSize',14);
plot(-1.*[fe_log10a(1) fe_log10a(1)],log10([fe_td50(2) fe_td50(3)]),'k.-.','LineWidth',1);
plot(-1.*[fe_log10a(2) fe_log10a(3)],log10([fe_td50(1) fe_td50(1)]),'k.-.','LineWidth',1);
h_re_meta=plot(-re_log10a(1),log10(re_td50(1)),'mx','LineWidth',2,'MarkerSize',14);
plot(-1.*[re_log10a(1) re_log10a(1)],log10([re_td50(2) re_td50(3)]),'m.-.','LineWidth',1);
plot(-1.*[re_log10a(2) re_log10a(3)],log10([re_td50(1) re_td50(1)]),'m.-.','LineWidth',1);

lgnd = legend([h_comb h_fe_meta h_re_meta],...
    'Pooled best fit',...
    ['Fixed-effects meta analysis',10,'best fit (projected)'],...
    ['Random-effects meta analysis',10,'best fit (projected)'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([-0.2 0.75]);
ylim([1.1 1.55]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('log_{10}(n)','FontSize',22); 
ylabel('log(TD50) (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_log10_n_log10_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_log10_n_log10_td50.png...']);
end


llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,comb_m_idx,:));

cur_fig=figure(50); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(1./(CGcomb.mLymanN),CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(1./(CGcomb.mLymanN(comb_n_idx)),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_a(1),fe_td50(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_a(1) fe_a(1)],[fe_td50(2) fe_td50(3)],'k.-.','LineWidth',1);
plot([fe_a(2) fe_a(3)],[fe_td50(1) fe_td50(1)],'k.-.','LineWidth',1);
h_re_meta=plot(re_a(1),re_td50(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_a(1) re_a(1)],[re_td50(2) re_td50(3)],'m.-.','LineWidth',1);
plot([re_a(2) re_a(3)],[re_td50(1) re_td50(1)],'m.-.','LineWidth',1);

lgnd = legend([h_comb h_fe_meta h_re_meta],...
    'Pooled best fit',...
    ['Fixed-effects meta analysis',10,'best fit (projected)'],...
    ['Random-effects meta analysis',10,'best fit (projected)'],...
    'Location','NorthWest');
set(lgnd,'FontSize',18);

xlim([0.1 1.5]);
ylim([12 34.5]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman a','FontSize',22); 
ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_a_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_a_td50.png...']);
end

%find position of meta n
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(comb_td50_idx,:,:));

cur_fig=figure(60); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(1./CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(1./CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_fe_meta=plot(fe_a(1),fe_m(1),'kx','LineWidth',2,'MarkerSize',14);
plot([fe_a(1) fe_a(1)],[fe_m(2) fe_m(3)],'k.-.','LineWidth',1);
plot([fe_a(2) fe_a(3)],[fe_m(1) fe_m(1)],'k.-.','LineWidth',1);
h_re_meta=plot(re_a(1),re_m(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_a(1) re_a(1)],[re_m(2) re_m(3)],'m.-.','LineWidth',1);
plot([re_a(2) re_a(3)],[re_m(1) re_m(1)],'m.-.','LineWidth',1);

lgnd = legend([h_comb h_fe_meta h_re_meta],...
    'Pooled best fit',...
    ['Fixed-effects meta analysis',10,'best fit (projected)'],...
    ['Random-effects meta analysis',10,'best fit (projected)'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([0.1 1.5]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman a','FontSize',22); 
ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_comb_a_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_comb_a_m.png...']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TD50 vs m 
% split at random meta analysis n
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,n_idx] = min(abs(CGcomb.mLymanN-re_n(1)));   
llhds = CGcomb.mLymanGrid.loglikelihood(:,:,n_idx);

cur_fig=figure(100); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(CGcomb.mLymanGrid.TD50,CGcomb.mLymanGrid.m,llhds',[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanGrid.TD50(comb_td50_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_re_meta=plot(re_td50(1),re_m(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_td50(1) re_td50(1)],[re_m(2) re_m(3)],'m.-','LineWidth',1);
plot([re_td50(2) re_td50(3)],[re_m(1) re_m(1)],'m.-','LineWidth',1);

lgnd = legend([h_comb h_re_meta],'Pooled best fit (projected)',['Random-effects',10,'Meta analysis best fit'],...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([15 32]);
ylim([0.04 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('TD_{50} (Gy)','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_re_meta_td50_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_re_meta_td50_m.png...']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot TD50 vs a 
% split at meta analysis m
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-re_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(200); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(CGcomb.mLymanN,CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',14);
h_re_meta=plot(re_n(1),re_td50(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_n(1) re_n(1)],[re_td50(2) re_td50(3)],'m.-','LineWidth',1);
plot([re_n(2) re_n(3)],[re_td50(1) re_td50(1)],'m.-','LineWidth',1);

lgnd = legend([h_comb h_re_meta],...
    'Pooled best fit (projected)',...
    ['Random-effects',10,'meta analysis best fit'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([0.26 3.2]);
ylim([12 34.5]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_re_meta_n_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_re_meta_n_td50.png...']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot m vs a 
% split at meta analysis TD50
%%%%%%%%%%%%%%%%%%%%%%%%%%

%find position of meta n
[~,td50_idx] = min(abs(CGcomb.mLymanGrid.TD50-re_td50(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(td50_idx,:,:));

cur_fig=figure(300); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_re_meta=plot(re_n(1),re_m(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_n(1) re_n(1)],[re_m(2) re_m(3)],'m.-','LineWidth',1);
plot([re_n(2) re_n(3)],[re_m(1) re_m(1)],'m.-','LineWidth',1);

lgnd = legend([h_comb h_re_meta],...
    'Pooled best fit (projected)',...
    ['Random-effects',10,'meta analysis best fit'],...
    'Location','SouthEast');
set(lgnd,'FontSize',18);

xlim([0.26 4.4]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman n','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_re_meta_n_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_re_meta_n_m.png...']);
end

%find position of meta n
[~,m_idx] = min(abs(CGcomb.mLymanGrid.m-re_m(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(:,m_idx,:));

cur_fig=figure(400); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(log10(CGcomb.mLymanN),log10(CGcomb.mLymanGrid.TD50),llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(log10(CGcomb.mLymanN(comb_n_idx)),log10(CGcomb.mLymanGrid.TD50(comb_td50_idx)),'r+','LineWidth',2,'MarkerSize',14);
h_re_meta=plot(-re_log10a(1),log10(re_td50(1)),'mx','LineWidth',2,'MarkerSize',14);
plot(-1.*[re_log10a(1) re_log10a(1)],log10([re_td50(2) re_td50(3)]),'m.-','LineWidth',1);
plot(-1.*[re_log10a(2) re_log10a(3)],log10([re_td50(1) re_td50(1)]),'m.-','LineWidth',1);

lgnd = legend([h_comb h_re_meta],...
    'Pooled best fit (projected)',...
    ['Random-effects',10,'meta analysis best fit'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([-0.2 0.75]);
ylim([1.1 1.55]);

set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('log_{10}(n)','FontSize',22); ylabel('log(TD50) (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_re_meta_log10_n_log10_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_re_meta_log10_n_log10_td50.png...']);
end


cur_fig=figure(500); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(1./(CGcomb.mLymanN),CGcomb.mLymanGrid.TD50,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(1./(CGcomb.mLymanN(comb_n_idx)),CGcomb.mLymanGrid.TD50(comb_td50_idx),'r+','LineWidth',2,'MarkerSize',14);
h_re_meta=plot(re_a(1),re_td50(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_a(1) re_a(1)],[re_td50(2) re_td50(3)],'m.-','LineWidth',1);
plot([re_a(2) re_a(3)],[re_td50(1) re_td50(1)],'m.-','LineWidth',1);

lgnd = legend([h_comb h_re_meta],...
    'Pooled best fit (projected)',...
    ['Random-effects',10,'meta analysis best fit'],...
    'Location','NorthWest');
set(lgnd,'FontSize',18);

xlim([0.1 1.5]);
ylim([12 34.5]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman a','FontSize',22); 
ylabel('TD50 (Gy)','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_re_meta_a_td50'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_re_meta_a_td50.png...']);
end


%find position of meta n
[~,td50_idx] = min(abs(CGcomb.mLymanGrid.TD50-re_td50(1)));   
llhds = squeeze(CGcomb.mLymanGrid.loglikelihood(td50_idx,:,:));

cur_fig=figure(600); clf reset;
set(cur_fig,'Position',ss_four2three);

[~,Htmp]=contour(1./CGcomb.mLymanN,CGcomb.mLymanGrid.m,llhds,[low99,low95,low68],'LineWidth',2);hold on; 
Cld = get(Htmp,'Children');
for j=1:length(Cld)
    if strcmp(get(Cld(j),'Type'),'patch')
        set(Cld(j),'EdgeColor',cntr_colors{j+(3-length(Cld))})
    end
end
h_comb=plot(1./CGcomb.mLymanN(comb_n_idx),CGcomb.mLymanGrid.m(comb_m_idx),'r+','LineWidth',2,'MarkerSize',14);
h_re_meta=plot(re_a(1),re_m(1),'mx','LineWidth',2,'MarkerSize',14);
plot([re_a(1) re_a(1)],[re_m(2) re_m(3)],'m.-','LineWidth',1);
plot([re_a(2) re_a(3)],[re_m(1) re_m(1)],'m.-','LineWidth',1);

lgnd = legend([h_comb h_re_meta],...
    'Pooled best fit (projected)',...
    ['Random-effects',10,'meta analysis best fit'],...
    'Location','NorthEast');
set(lgnd,'FontSize',18);

xlim([0.1 1.5]);
ylim([0.12 0.65]);
 
set(gca,'xminortick','on','yminortick','on');
set(gca,'box','on');
set(gca,'FontSize',18)
xlabel('Lyman a','FontSize',22); ylabel('Lyman m','FontSize',22);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'map_llhd_re_meta_a_m'],'-png');
    disp(['Saving ',fig_loc,'map_llhd_re_meta_a_m.png...']);
end
end