function EudAnalysis_MetaAnalysis_Compare
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
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';

load(fn,'CGcomb');

% change to log10(n)  (add - for log10(a)

LymanN = log10(CGcomb.mLymanN);
CGcomb.mLymanN = LymanN;
td50s = CGcomb.mLymanGridTD50Range;
ms = CGcomb.mLymanGridMRange;
% Want to find likelihoods for meta-analysis parameters using pooled
% parameter grid

% llhds (TD50,m,n)
llhds = CGcomb.mLymanGrid.loglikelihood;
clear CGcomb;

%find TD50 ind
[~,td50_ind] = min(abs(td50s-re_td50(1)));
[~,td50_ll_ind] = min(abs(td50s-re_td50(2)));
[~,td50_ul_ind] = min(abs(td50s-re_td50(3)));
%find m_ind 
[~,m_ind] = min(abs(ms-re_m(1)));
[~,m_ll_ind] = min(abs(ms-re_m(2)));
[~,m_ul_ind] = min(abs(ms-re_m(3)));
%find a_ind
[~,a_ind] = min(abs(-LymanN-re_log10a(1)));
[~,a_ll_ind] = min(abs(-LymanN-re_log10a(2)));
[~,a_ul_ind] = min(abs(-LymanN-re_log10a(3)));

re_llhd = llhds(td50_ind,m_ind,a_ind);
re_ll_llhd = llhds(td50_ll_ind,m_ll_ind,a_ll_ind);
re_ul_llhd = llhds(td50_ul_ind,m_ul_ind,a_ul_ind);

%% pooled

[mx_llhd,loc] = max(llhds(:));
[pld_td50_ind,pld_m_ind,pld_n_ind] = ind2sub(size(llhds),loc);

%% L >= L(max) - chi2inv(alpha,n(parameters)/2;
% alpha = chi2cdf(-2*(L - L(max)),n(parameters));
pval = 1-chi2cdf(-2*(re_llhd-mx_llhd),3);
pval_ll = 1-chi2cdf(-2*(re_ll_llhd-mx_llhd),3);
pval_ul = 1-chi2cdf(-2*(re_ul_llhd-mx_llhd),3);

disp(['Pooled gEUD log-likelihood: ',num2str(llhds(pld_td50_ind,pld_m_ind,pld_n_ind)),' test: ',num2str(mx_llhd)]);
disp(['Random-Effects log-likelihood: ',num2str(re_llhd)]);

disp([]);
disp(['Random-Effects exclusion p-value is: ',num2str(pval,3)]);



end