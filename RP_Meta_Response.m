function RP_Meta_Response

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
fig_loc = 'Z:\elw\MATLAB\meta_analy\slides\figures\meta_response\latest\';

do_print = true;
cur_fig_ctr = 1;


% pld_vals = [0.69,1.45,19.2,0.38];
% pld_lcl = [0.33,0.94,14.3,0.29];
% pld_ucl = [1.06,3.03,26.2,0.51];

% loglikelihoods extracted from analysis objet in RP_Write_LLHDs.m
fn='Z:/elw/MATLAB/meta_analy/meta_data/lkb_comb_llhds.mat';
load(fn);% lkb_llhds, lkb_n, lkb_m, lkb_td50


%% from one-sided meta
% re_a = [1.095,0.4141,1.776];
% % re_n = [1.071,0.1658,1.976]
% re_td50 = [24.24,14.46,34.03];
% re_m = [0.3438,0.1412,0.5463];
% 
% fe_a = [1.103,0.6382,1.567];
% % fe_n = [0.8603,0.4404,1.28]
% fe_m = [0.3066,0.2221,0.3911];
% fe_td50 = [22.94,17.24,28.63];

%% from two-sided meta
fe_a = [0.8431,0.437,1.249];
% fe_n = [0.757,0.4388,1.075];
fe_td50 = [21.89,15.96,27.81] ;
fe_m = [0.2144,0.14,0.2888];

re_a = [0.8035,0.1965,1.411];
% re_n = [1.239,0.04706,2.431];
re_td50 = [21.89,15.96,27.81];
re_m = [0.287,0.1311,0.4428];

% a=0.69;
cur_re_a = re_a(1);
[~,re_a_ind] = min(abs((1./lkb_n)-cur_re_a)); 
ll = lkb_llhds(:,:,re_a_ind); % log likelihood of log10(a) = loga
re_mx = max(ll(:));
[re_xx,re_yy] = find(ll == re_mx); % the coefficients
cur_re_td50 = lkb_td50(re_xx);
cur_re_m = lkb_m(re_yy);
disp(['Best Pooled RE (td50,m,a) = (',num2str([cur_re_td50,cur_re_m,cur_re_a]),')']);
pld_re_vals = [cur_re_a, 1/cur_re_a, cur_re_td50, cur_re_m];


cur_fe_a = fe_a(1);
[~,fe_a_ind] = min(abs((1./lkb_n)-cur_fe_a)); 
ll = lkb_llhds(:,:,fe_a_ind); % log likelihood of log10(a) = loga
fe_mx = max(ll(:));
[fe_xx,fe_yy] = find(ll == fe_mx); % the coefficients
cur_fe_td50 = lkb_td50(fe_xx);
cur_fe_m = lkb_m(fe_yy);
disp(['Best Pooled FE (td50,m,a) = (',num2str([cur_fe_td50,cur_fe_m,cur_fe_a]),')']);
pld_fe_vals = [cur_fe_a, 1/cur_fe_a, cur_fe_td50, cur_fe_m];




doses = [0:0.1:25]';

%% pooled
% lyman probability
lkb_re_prob = normcdf((doses-pld_re_vals(3))./(pld_re_vals(3)*pld_re_vals(4)),0,1);
% low68 = mx - 0.5*(1*2);
lkb_re_llhd_low68 = re_mx - 0.5*(chi2inv(0.68,3));
lkb_re_low68_td50_m_ctr = contourc(lkb_td50,lkb_m,ll',[lkb_re_llhd_low68 lkb_re_llhd_low68]);
lkb_re_low68_td50 = lkb_re_low68_td50_m_ctr(1,2:end);
lkb_re_low68_m = lkb_re_low68_td50_m_ctr(2,2:end);

lkb_re_td50_68_ll = min(lkb_re_low68_td50); 
lkb_re_td50_68_ul = max(lkb_re_low68_td50); 

lkb_re_m_68_ll = min(lkb_re_low68_m); 
lkb_re_m_68_ul = max(lkb_re_low68_m); 

% Random-Effects
re_prob = normcdf((doses-re_td50(1))./(re_m(1).*re_td50(1)),0,1);

re_td50_lsd = (re_td50(1)-re_td50(2))./1.96;
re_td50_usd = (re_td50(3)-re_td50(1))./1.96;
re_td50_dist_ll = normrnd(re_td50(1),re_td50_lsd,1e5,1);
re_td50_dist_ul = normrnd(re_td50(1),re_td50_usd,1e5,1);
re_td50_dist_68_ll = re_td50(1)-re_td50_lsd;
re_td50_dist_68_ul = re_td50(1)+re_td50_usd;


re_m_lsd = (re_m(1)-re_m(2))./1.96;
re_m_usd = (re_m(3)-re_m(1))./1.96;
re_m_dist_ll = normrnd(re_m(1),re_m_lsd,1e5,1);
re_m_dist_ul = normrnd(re_m(1),re_m_usd,1e5,1);
re_m_dist_68_ll = re_m(1)-re_m_lsd;
re_m_dist_68_ul = re_m(1)+re_m_usd;

% Fixed-Effects
fe_prob = normcdf((doses-fe_td50(1))./(fe_m(1).*fe_td50(1)),0,1);

fe_td50_lsd = (fe_td50(1)-fe_td50(2))./1.96;
fe_td50_usd = (fe_td50(3)-fe_td50(1))./1.96;
fe_td50_dist_ll = normrnd(fe_td50(1),fe_td50_lsd,1e5,1);
fe_td50_dist_ul = normrnd(fe_td50(1),fe_td50_usd,1e5,1);
fe_td50_dist_68_ll = fe_td50(1)-fe_td50_lsd;
fe_td50_dist_68_ul = fe_td50(1)+fe_td50_usd;

fe_m_lsd = (fe_m(1)-fe_m(2))./1.96;
fe_m_usd = (fe_m(3)-fe_m(1))./1.96;
fe_m_dist_ll = normrnd(fe_m(1),fe_m_lsd,1e5,1);
fe_m_dist_ul = normrnd(fe_m(1),fe_m_usd,1e5,1);
fe_m_dist_68_ll = fe_m(1)-fe_m_lsd;
fe_m_dist_68_ul = fe_m(1)+fe_m_usd;

% FE lyman probability
lkb_fe_prob = normcdf((doses-pld_fe_vals(3))./(pld_fe_vals(3)*pld_fe_vals(4)),0,1);
% low68 = mx - 0.5*(1*2);
lkb_fe_llhd_low68 = fe_mx - 0.5*(chi2inv(0.68,3));
lkb_fe_low68_td50_m_ctr = contourc(lkb_td50,lkb_m,ll',[lkb_fe_llhd_low68 lkb_fe_llhd_low68]);
lkb_fe_low68_td50 = lkb_fe_low68_td50_m_ctr(1,2:end);
lkb_fe_low68_m = lkb_fe_low68_td50_m_ctr(2,2:end);

lkb_fe_td50_68_ll = min(lkb_fe_low68_td50); 
lkb_fe_td50_68_ul = max(lkb_fe_low68_td50); 

lkb_fe_m_68_ll = min(lkb_fe_low68_m); 
lkb_fe_m_68_ul = max(lkb_fe_low68_m); 


% at each does point, sample parameter CIs, pick largest fluctuation
re_prob_68_lcl = inf(length(doses),1);
re_prob_68_ucl = inf(length(doses),1);

fe_prob_68_lcl = inf(length(doses),1);
fe_prob_68_ucl = inf(length(doses),1);

lkb_re_prob_68_lcl = zeros(length(doses),1);
lkb_re_prob_68_ucl = zeros(length(doses),1);

lkb_fe_prob_68_lcl = zeros(length(doses),1);
lkb_fe_prob_68_ucl = zeros(length(doses),1);

for d=1:length(doses)
    % RE pooled
    lkb_re_prob_68 = normcdf((doses(d)-lkb_re_low68_td50)./(lkb_re_low68_m.*lkb_re_low68_td50),0,1); % probability at each iso point
    lkb_re_prob_68_lcl(d) = min(lkb_re_prob_68);
    lkb_re_prob_68_ucl(d) = max(lkb_re_prob_68);
    
    
    % Random-effects
    re_prob_ul = normcdf((doses(d)-re_td50_dist_ul)./(re_m_dist_ul.*re_td50_dist_ul),0,1);
    re_std_ul = std(re_prob_ul);
    re_prob_ll = normcdf((doses(d)-re_td50_dist_ll)./(re_m_dist_ll.*re_td50_dist_ll),0,1);
    re_std_ll = std(re_prob_ll);
  
    re_prob_68_lcl(d) = max((re_prob(d)-(re_std_ll)),0);
    re_prob_68_ucl(d) = min((re_prob(d)+(re_std_ul)),1);
    
    % Fixed-effects
    fe_prob_ul = normcdf((doses(d)-fe_td50_dist_ul)./(fe_m_dist_ul.*fe_td50_dist_ul),0,1);
    fe_std_ul = std(fe_prob_ul);
    fe_prob_ll = normcdf((doses(d)-fe_td50_dist_ll)./(fe_m_dist_ll.*fe_td50_dist_ll),0,1);
    fe_std_ll = std(fe_prob_ll);
  
    fe_prob_68_lcl(d) = max((fe_prob(d)-(fe_std_ll)),0);
    fe_prob_68_ucl(d) = min((fe_prob(d)+(fe_std_ul)),1);
    
    % FE pooled
    lkb_fe_prob_68 = normcdf((doses(d)-lkb_fe_low68_td50)./(lkb_fe_low68_m.*lkb_fe_low68_td50),0,1); % probability at each iso point
    lkb_fe_prob_68_lcl(d) = min(lkb_fe_prob_68);
    lkb_fe_prob_68_ucl(d) = max(lkb_fe_prob_68);
        
    
end



cur_fig=figure(cur_fig_ctr); clf reset; hold on;
set(cur_fig,'Position',ss_four2three);
cur_fig_ctr = cur_fig_ctr+1;

h_lkb_re_prob=plot(doses,lkb_re_prob,'r-','LineWidth',2);hold on;
plot(doses,lkb_re_prob_68_ucl,'r-.','LineWidth',1);
plot(doses,lkb_re_prob_68_lcl,'r-.','LineWidth',1);
h_re_prob=plot(doses,re_prob,'k-','LineWidth',2);
plot(doses,re_prob_68_ucl,'k-.','LineWidth',1);
plot(doses,re_prob_68_lcl,'k-.','LineWidth',1);

lgnd = legend([h_re_prob h_lkb_re_prob],...
    ['\underline{Random-Effect} [68\% CIs]',10,...
    '~~TD$_{50} = ',num2str(re_td50(1),3),...
    '~[',num2str(re_td50_dist_68_ll,3),'-',num2str(re_td50_dist_68_ul,3),']$',10,...
    '~~$m = ',num2str(re_m(1),3),...
    '~[',num2str(re_m_dist_68_ll,3),'-',num2str(re_m_dist_68_ul,3),']$'],...
    ['\underline{Pooled} [68\% CIs]',10,...
    '~~TD$_{50} = ',num2str(pld_re_vals(3),3),...
                  '~[',num2str(lkb_re_td50_68_ll,3),'-',num2str(lkb_re_td50_68_ul,3),']$',10,...
                  '~~$m = ',num2str(pld_re_vals(4),3),...
                  '~[',num2str(lkb_re_m_68_ll,3),'-',num2str(lkb_re_m_68_ul,3),']$'],...
    'location','NorthWest');

set(lgnd,'interpreter','latex');            
set(lgnd,'FontSize',22);

% text(0.05,0.7,['\underline{Random-Effect}:',10,...
%                   '~~TD$_{50} = ',num2str(re_td50(1),3),...
%                   '~[',num2str(re_td50_dist_68_ll,3),'-',num2str(re_td50_dist_68_ul,3),']$',10,...
%                   '~~$m = ',num2str(re_m(1),3),...
%                   '~[',num2str(re_m_dist_68_ll,3),'-',num2str(re_m_dist_68_ul,3),']$',10,...
%                    '\underline{Pooled}:',10,...
%                   '~~TD$_{50} = ',num2str(pld_re_vals(3),3),...
%                   '~[',num2str(lkb_re_td50_68_ll,3),'-',num2str(lkb_re_td50_68_ul,3),']$',10,...
%                   '~~$m = ',num2str(pld_re_vals(4),3),...
%                   '~[',num2str(lkb_re_m_68_ll,3),'-',num2str(lkb_re_m_68_ul,3),']$'],...
%                   'unit','normalize',...
%               'interpreter','latex',...
%                 'fontsize',18);

% '{\color[rgb]{1 0 0}Pooled}:',10,...            
            
set(gca,'xminortick','on','yminortick','on');

set(gca,'FontSize',18);
xlabel(['gEUD (a=',num2str(pld_re_vals(1),3),') [Gy]'],'FontSize',20);
ylabel('RP probability','FontSize',20);
ylim([0 0.8]);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'lkb_re_response'],'-pdf');
    disp(['Saving ',fig_loc,'lkb_re_response.pdf...']);
end

cur_fig=figure(cur_fig_ctr); clf reset; hold on;
set(cur_fig,'Position',ss_four2three);
cur_fig_ctr = cur_fig_ctr+1;

h_lkb_fe_prob=plot(doses,lkb_fe_prob,'r-','LineWidth',2);hold on;
plot(doses,lkb_fe_prob_68_ucl,'r-.','LineWidth',1);
plot(doses,lkb_fe_prob_68_lcl,'r-.','LineWidth',1);
h_fe_prob=plot(doses,fe_prob,'k-','LineWidth',2);
plot(doses,fe_prob_68_ucl,'k-.','LineWidth',1);
plot(doses,fe_prob_68_lcl,'k-.','LineWidth',1);

lgnd = legend([h_fe_prob h_lkb_fe_prob],...
    ['\underline{Fixed-Effect} [68\% CIs]',10,...
    '~~TD$_{50} = ',num2str(fe_td50(1),3),...
    '~[',num2str(fe_td50_dist_68_ll,3),'-',num2str(fe_td50_dist_68_ul,3),']$',10,...
    '~~$m = ',num2str(fe_m(1),3),...
    '~[',num2str(fe_m_dist_68_ll,3),'-',num2str(fe_m_dist_68_ul,3),']$'],...
    ['\underline{Pooled} [68\% CIs]',10,...
    '~~TD$_{50} = ',num2str(pld_fe_vals(3),3),...
                  '~[',num2str(lkb_fe_td50_68_ll,3),'-',num2str(lkb_fe_td50_68_ul,3),']$',10,...
                  '~~$m = ',num2str(pld_fe_vals(4),3),...
                  '~[',num2str(lkb_fe_m_68_ll,3),'-',num2str(lkb_fe_m_68_ul,3),']$'],...
    'location','NorthWest');

set(lgnd,'interpreter','latex');            
set(lgnd,'FontSize',22);
% 
% text(0.05,0.7,['\underline{Fixed-Effect}:',10,...
%                   '~~TD$_{50} = ',num2str(fe_td50(1),3),...
%                   '~[',num2str(fe_td50_dist_68_ll,3),'-',num2str(fe_td50_dist_68_ul,3),']$',10,...
%                   '~~$m = ',num2str(fe_m(1),3),...
%                   '~[',num2str(fe_m_dist_68_ll,3),'-',num2str(fe_m_dist_68_ul,3),']$',10,...
%                    '\underline{Pooled}:',10,...
%                   '~~TD$_{50} = ',num2str(pld_fe_vals(3),3),...
%                   '~[',num2str(lkb_fe_td50_68_ll,3),'-',num2str(lkb_fe_td50_68_ul,3),']$',10,...
%                   '~~$m = ',num2str(pld_fe_vals(4),3),...
%                   '~[',num2str(lkb_fe_m_68_ll,3),'-',num2str(lkb_fe_m_68_ul,3),']$'],...
%                   'unit','normalize',...
%               'interpreter','latex',...
%                 'fontsize',18);

set(gca,'xminortick','on','yminortick','on');

set(gca,'FontSize',18);
xlabel(['gEUD (a=',num2str(pld_fe_vals(1),5),') [Gy]'],'FontSize',20);
ylabel('RP probability','FontSize',20);
ylim([0 0.8]);

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'lkb_fe_response'],'-pdf');
    disp(['Saving ',fig_loc,'lkb_fe_response.pdf...']);
end


end