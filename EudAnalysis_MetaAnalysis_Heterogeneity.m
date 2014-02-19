function EudAnalysis_MetaAnalysis_Heterogeneity
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