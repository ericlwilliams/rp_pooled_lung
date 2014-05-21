function RP_Write_LLHDs
    % writes file of log-likelihood matrix only
% used in RP_InterInst_Response (eg) for CI calculations
    fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
    load(fn,'CGcomb');
    
    lkb_llhds = CGcomb.mLymanGrid.loglikelihood;
    lkb_n = CGcomb.mLymanN;
    lkb_td50 = CGcomb.mLymanGrid.TD50;
    lkb_m = CGcomb.mLymanGrid.m;

    fname='Z:/elw/MATLAB/meta_analy/meta_data/lkb_comb_llhds.mat';
    save(fname,'lkb_llhds','lkb_n','lkb_td50','lkb_m');
    
    
end