function EudAnalysis_MetaAnalysis_Heterogeneity
% template from EudAnalysis_MetaAnalysis_Compare
tic; %close all;

%% random effects meta analysis results [val,95% lo, hi]
%re_td50 = [21.88,15.9,27.8];
%re_log10a = [-0.0862,-0.69897,0.15533];
%re_m = [0.29,0.13,0.44];

%% fixed effects (even though labelled as re_)
re_td50 = [21.88,15.9,27.8];
re_log10a = [-0.0506,-0.2924,0.1038];
re_m = [0.21,0.14,0.28];

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';

%% COMB only
 fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
% fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta.mat';
%fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_crs_EUD_fine_meta.mat';

if isunix %on mac
    fig_loc = strrep(fig_loc,'Z:/elw/','/Users/elw/Documents/');
    fn = strrep(fn,'Z:/elw/','/Users/elw/Documents/');
end

load(fn,'CGnki','CGmsk','CGrtog','CGum','CGcomb');
protocols = {'NKI','MSK','RTOG','UM','Comb'};
% change to log10(n)  (add - for log10(a)
% LymanN = log10(CGmsk.mLymanN);
% CGcomb.mLymanN = LymanN;
% td50s = CGcomb.mLymanGridTD50Range;
% ms = CGcomb.mLymanGridMRange;
% Want to find likelihoods for meta-analysis parameters using pooled
% parameter grid

% llhds (TD50,m,n)
% nki_llhds = CGnki.mLymanGrid.loglikelihood;
% rtog_llhds = CGrtog.mLymanGrid.loglikelihood;
% umich_llhds = CGum.mLymanGrid.loglikelihood;
% msk_llhds = CGmsk.mLymanGrid.loglikelihood;
% clear CGcomb;clear CGnki;clear CGmsk;clear CGrtog;

% Quantifying heterogeneity with modified (?) chi2 test using likelihoods
% Need mean values of parameters first
  
% [mx_loga,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
%     [~,~,a_loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
%     loga = -OCobj.mLymanN(a_loc);
%     
    
OCobjs = [CGnki;CGmsk;CGrtog;CGum;CGcomb]';


% prepare
num = length(OCobjs(:));
% Q = zeros(3,1); Qp = zeros(3,1);
% I2 = zeros(3,1); I2up = zeros(3,1); I2down = zeros(3,1);
b0 = zeros(num,1); b1 = zeros(num,1); a = zeros(num,1);
b0sd = zeros(num,1); b1sd = zeros(num,1); asd = zeros(num,1);

% best parameters and standard deviations
% 68% CI found for each parameter, SD taken from those.
% TODO, use CIs to get I2 range?
for k = 1:num
    % best parameters
    mx = max(OCobjs(k).mLymanGrid.loglikelihood(:));
    
    %OCobj.mLymanGrid.TD50

    % HERE
    likely1 = max(max(OCobjs(k).mLymanGrid.loglikelihood,[],3),[],2);
    likely2 = max(max(OCobjs(k).mLymanGrid.loglikelihood,[],3),[],1);
    likely3 = max(max(OCobjs(k).mLymanGrid.loglikelihood,[],2),[],1);
    
    %CI68 = mx - 0.5*1; % 68% confidence interval formula, *1
    %because 68% is 1 sigma away from the peak
    %tmp
    CI68 = mx - 0.5*1; % 68% confidence interval formula, *1 because 68% is 1 sigma away from the peak
    % 68% point of b0
    [b0(k), CI] = ConfidenceInterval(OCobjs(k).mLymanGrid.TD50,likely1, CI68);
    b0sd(k) = min(CI(2)-b0(k), b0(k)-CI(1));
    % 68% point of b1
    [b1(k), CI] = ConfidenceInterval(OCobjs(k).mLymanGrid.m,likely2, CI68);
    b1sd(k) = min(CI(2)-b1(k), b1(k)-CI(1));
    % 68% point of a
    [a(k), CI] = ConfidenceInterval(OCobjs(k).mLymanN,likely3, CI68);
    asd(k) = min(CI(2)-a(k), a(k)-CI(1));
   
end
   
% Get mean values of parameters, likelihood ratio between mean value and 
% meta_n = mean(a);
% meta_td50 = mean(b0);
% meta_m = mean(b1);


meta_td50 = 21.88;% random effects
meta_n = 1.22;% random effects
meta_m = 0.29; % random effects


disp(['Meta: ',num2str(meta_td50),',',num2str(meta_m),',',num2str(meta_n)]);
pseudo_pvals = zeros(num,1);
for j = 1:num
    mx = max(OCobjs(j).mLymanGrid.loglikelihood(:));
    % for each inst. find likelihood at meta values
    [~,td50_ind] = min(abs(OCobjs(j).mLymanGrid.TD50 - meta_td50));
    [~,m_ind] = min(abs(OCobjs(j).mLymanGrid.m - meta_m));
    [~,n_ind] = min(abs(OCobjs(j).mLymanN - meta_n));
    meta_llhd = OCobjs(j).mLymanGrid.loglikelihood(td50_ind,m_ind,n_ind);
    
    pseudo_chi2stat = -2*meta_llhd+2*mx;
    pseudo_pvals(j) = 1-chi2cdf(pseudo_chi2stat,3);
    
    disp([protocols{j},'(',num2str(b0(j)),',',num2str(b1(j)),',',num2str(a(j)),')']);
    disp([protocols{j},' LLHD: ',num2str(mx)]);
    disp(['Meta LLHD: ',num2str(meta_llhd)]);
    disp(['Chi2: ',num2str(pseudo_chi2stat)]);
    disp(['P-val: ',num2str(pseudo_pvals(j))]);
    
end
   

disp(pseudo_pvals)
    

%% Compare protocol results from combined llhd matrix

disp(['Compare protocol results from combined llhd matrix' ]);
comb_pseudo_pvals = zeros(num-1,1);
mx = max(OCobjs(5).mLymanGrid.loglikelihood(:));
for j = 1:(num-1)
    cur_td50 = b0(j);
    cur_m = b1(j);
    cur_n = a(j);
    % for each inst. find likelihood at meta values
    [~,td50_ind] = min(abs(OCobjs(5).mLymanGrid.TD50 - cur_td50));
    [~,m_ind] = min(abs(OCobjs(5).mLymanGrid.m - cur_m));
    [~,n_ind] = min(abs(OCobjs(5).mLymanN - cur_n));
    meta_llhd = OCobjs(5).mLymanGrid.loglikelihood(td50_ind,m_ind,n_ind);
    
    pseudo_chi2stat = -2*meta_llhd+2*mx;
    comb_pseudo_pvals(j) = 1-chi2cdf(pseudo_chi2stat,3);
    
    disp([protocols{j},'(',num2str(b0(j)),',',num2str(b1(j)),',',num2str(a(j)),')']);
    disp([protocols{j},' LLHD: ',num2str(meta_llhd)]);
    disp(['Meta LLHD: ',num2str(mx)]);
    disp(['Chi2: ',num2str(pseudo_chi2stat)]);
    disp(['P-val: ',num2str(comb_pseudo_pvals(j))]);
    
end

% meta analysis model in pooled matrix
  [~,td50_ind] = min(abs(OCobjs(5).mLymanGrid.TD50 - meta_td50));
    [~,m_ind] = min(abs(OCobjs(5).mLymanGrid.m - meta_m));
    [~,n_ind] = min(abs(OCobjs(5).mLymanN - meta_n));
    meta_llhd = OCobjs(5).mLymanGrid.loglikelihood(td50_ind,m_ind,n_ind);
    
    pseudo_chi2stat = -2*meta_llhd+2*mx;
    meta_pseudo_pvals = 1-chi2cdf(pseudo_chi2stat,3);


disp(comb_pseudo_pvals)

disp(['meta paramters in pooled matrix p-value: ',num2str(meta_pseudo_pvals)]);


