function RP_InterInst_Response
% Evaluate parameter uncertanties by incorporating inter-institutional
% variances (from rpMetaAnalysis)

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
fig_loc = 'Z:\elw\MATLAB\meta_analy\slides\figures\lkb_response\latest\';

do_print = true;
cur_fig_ctr = 1;


%% don't use full parameter 95% CIs...
% we are setting a = 0.69, need to find 
% CIs only in that plane!

%   a,    n, td50,  m
% pld_vals = [0.69,1.45,19.2,0.38];
% pld_lcl = [0.33,0.94,14.3,0.29];
% pld_ucl = [1.06,3.03,26.2,0.51];

% loglikelihoods extracted from analysis objet in RP_Write_LLHDs.m
fn='Z:/elw/MATLAB/meta_analy/meta_data/lkb_comb_llhds.mat';
load(fn);% lkb_llhds, lkb_n, lkb_m, lkb_td50

ii_t2_m = 0.01;

a=0.69;
[~,a_ind] = min(abs((1./lkb_n)-a)); 
ll = lkb_llhds(:,:,a_ind); % log likelihood of log10(a) = loga
mx = max(ll(:));
[xx,yy] = find(ll == mx); % the coefficients
td50 = lkb_td50(xx);
m = lkb_m(yy);
disp(['Best (td50,m,a) = (',num2str([td50,m,a]),')']);
pld_vals = [a, 1/a, td50, m];

doses = [0:0.1:25]';
td50_step = 0.1;

%% 68% CI 
% low68 = mx - 0.5*(1*2);
low68 = mx - 0.5*(chi2inv(0.68,2));
low68_td50_m_ctr = contourc(lkb_td50,lkb_m,ll',[low68 low68]);
low68_td50 = low68_td50_m_ctr(1,2:end);
low68_m = low68_td50_m_ctr(2,2:end);

% step through td50, find CI on m, calc sigma for m, add t2.
low68_td50_bins = [min(low68_td50)+td50_step:td50_step:max(low68_td50)+td50_step];

low68_ii_m_min = inf(1,length(low68_td50_bins));
low68_ii_m_max = inf(1,length(low68_td50_bins));

low68_m_min = inf(1,length(low68_td50_bins));
low68_m_max = inf(1,length(low68_td50_bins));

low68_m_sig = inf(1,length(low68_td50_bins));
low68_ii_m_sig = inf(1,length(low68_td50_bins));


for i=1:length(low68_td50_bins)
   
    cur_low68_inds = [low68_td50<low68_td50_bins(i) & low68_td50>=(low68_td50_bins(i)-td50_step)];
    cur_low68_m = low68_m(cur_low68_inds);
    
    low68_m_min(i) = min(cur_low68_m);
    low68_m_max(i) = max(cur_low68_m);

    low68_m_sig(i) = (low68_m_max(i)-low68_m_min(i))/2;
    low68_ii_m_sig(i) = sqrt(low68_m_sig(i)^2 + ii_t2_m);
    
    low68_ii_m_max(i) = median(cur_low68_m) + low68_ii_m_sig(i);
    low68_ii_m_min(i) = median(cur_low68_m) - low68_ii_m_sig(i);
    
end
% re-center
low68_td50_bins = low68_td50_bins-(td50_step/2);

%% 95% CI 
% low95 = mx - 0.5*(1*2);
low95 = mx - 0.5*(chi2inv(0.95,2));
low95_td50_m_ctr = contourc(lkb_td50,lkb_m,ll',[low95 low95]);
low95_td50 = low95_td50_m_ctr(1,2:end);
low95_m = low95_td50_m_ctr(2,2:end);

low95_td50_bins = [min(low95_td50)+td50_step:td50_step:max(low95_td50)+td50_step];

low95_ii_m_min = inf(1,length(low95_td50_bins));
low95_ii_m_max = inf(1,length(low95_td50_bins));

low95_m_min = inf(1,length(low95_td50_bins));
low95_m_max = inf(1,length(low95_td50_bins));

for j=1:length(low95_td50_bins)
   
    cur_low95_inds = [low95_td50<low95_td50_bins(j) & low95_td50>=(low95_td50_bins(j)-td50_step)];
    cur_low95_m = low95_m(cur_low95_inds);
    
    low95_m_min(j) = min(cur_low95_m);
    low95_m_max(j) = max(cur_low95_m);

    cur_low95_m_sig = (low95_m_max(j)-low95_m_min(j))/(2*1.96);
    cur_low95_ii_m_sig = sqrt(cur_low95_m_sig^2 + ii_t2_m);
    
    low95_ii_m_max(j) = median(cur_low95_m) + (1.96*cur_low95_ii_m_sig);
    low95_ii_m_min(j) = median(cur_low95_m) - (1.96*cur_low95_ii_m_sig);
    
end
% re-center
low95_td50_bins = low95_td50_bins-(td50_step/2);





%% lyman probability
lkb_stat = (doses-pld_vals(3))./(pld_vals(4)*pld_vals(3));
lkb_prob = normcdf(lkb_stat,0,1);

lkb_prob_68_lcl = zeros(length(doses),1);
lkb_prob_68_ucl = zeros(length(doses),1);

lkb_prob_68_ii_lcl = zeros(length(doses),1);
lkb_prob_68_ii_ucl = zeros(length(doses),1);

lkb_prob_95_lcl = zeros(length(doses),1);
lkb_prob_95_ucl = zeros(length(doses),1);

lkb_prob_95_ii_lcl = zeros(length(doses),1);
lkb_prob_95_ii_ucl = zeros(length(doses),1);

for d=1:length(doses)
    
    lkb_prob_68_min = normcdf((doses(d)-low68_td50_bins)./(low68_m_min.*low68_td50_bins),0,1);
    lkb_prob_68_max = normcdf((doses(d)-low68_td50_bins)./(low68_m_max.*low68_td50_bins),0,1);
    
    lkb_prob_68_lcl(d) = max(min(min([lkb_prob_68_min;lkb_prob_68_max])),0);
    lkb_prob_68_ucl(d) = max(max([lkb_prob_68_min;lkb_prob_68_max]));
      
   % adding inter-inst heterogeneity
    lkb_prob_68_ii_min = normcdf((doses(d)-low68_td50_bins)./(low68_ii_m_min.*low68_td50_bins),0,1);
    lkb_prob_68_ii_max = normcdf((doses(d)-low68_td50_bins)./(low68_ii_m_max.*low68_td50_bins),0,1);
    
    lkb_prob_68_ii_lcl(d) = max(min(min([lkb_prob_68_ii_min;lkb_prob_68_ii_max])),0);
    lkb_prob_68_ii_ucl(d) = max(max([lkb_prob_68_ii_min;lkb_prob_68_ii_max]));

    
    
    lkb_prob_95_min = normcdf((doses(d)-low95_td50_bins)./(low95_m_min.*low95_td50_bins),0,1);
    lkb_prob_95_max = normcdf((doses(d)-low95_td50_bins)./(low95_m_max.*low95_td50_bins),0,1);
    
    lkb_prob_95_lcl(d) = max(min(min([lkb_prob_95_min;lkb_prob_95_max])),0);
    lkb_prob_95_ucl(d) = max(max([lkb_prob_95_min;lkb_prob_95_max]));
    
    
    lkb_prob_95_ii_min = normcdf((doses(d)-low95_td50_bins)./(low95_ii_m_min.*low95_td50_bins),0,1);
    lkb_prob_95_ii_max = normcdf((doses(d)-low95_td50_bins)./(low95_ii_m_max.*low95_td50_bins),0,1);
    
    lkb_prob_95_ii_lcl(d) = max(min(min([lkb_prob_95_ii_min;lkb_prob_95_ii_max])),0);
    lkb_prob_95_ii_ucl(d) = max(max([lkb_prob_95_ii_min;lkb_prob_95_ii_max]));

    
end

 
cur_fig=figure(cur_fig_ctr); clf reset; hold on;
set(cur_fig,'Position',ss_four2three);
cur_fig_ctr = cur_fig_ctr+1;

h_prob=plot(doses,lkb_prob,'k-','LineWidth',2);hold on;
h_68_prob=plot(doses,lkb_prob_68_ucl,'r-','LineWidth',1);
plot(doses,lkb_prob_68_lcl,'r-','LineWidth',1);
h_68_ii_prob=plot(doses,lkb_prob_68_ii_ucl,'r-.','LineWidth',2);
plot(doses,lkb_prob_68_ii_lcl,'r-.','LineWidth',2);
h_95_prob=plot(doses,lkb_prob_95_ucl,'b-','LineWidth',1);
plot(doses,lkb_prob_95_lcl,'b-','LineWidth',1);
h_95_ii_prob=plot(doses,lkb_prob_95_ii_ucl,'b-.','LineWidth',2);
plot(doses,lkb_prob_95_ii_lcl,'b-.','LineWidth',2);

lgnd=legend([h_68_prob h_68_ii_prob h_95_prob h_95_ii_prob],...
            '68\% CI','68\% CI ($\sigma^2 + \tau^2$)',...
            '95\% CI','95\% CI ($\sigma^2 + \tau^2$)',...
            'Location','NorthWest');
set(lgnd,'interpreter','latex');        
set(lgnd,'FontSize',18);

text(0.01,0.65,['Significant inter-institutional',10,...
        'heterogeneity in LKB',10,...
        'parameter $m$. From random-',10,...
        'effects meta-analysis,~$\tau^2_{m} = 0.01$'],...
        'interpreter','latex',...
        'unit','normalize',...
        'FontSize',18);
set(gca,'xminortick','on','yminortick','on');
% set(gca,'box','on');
set(gca,'FontSize',18);
xlabel('gEUD (a=0.69) [Gy]','FontSize',20);
ylabel('RP probability','FontSize',20);
ylim([0 0.8]);
if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'lkb_response'],'-pdf');
    disp(['Saving ',fig_loc,'lkb_response.pdf...']);
end


end