function  RP_Pooled_Results()
%RP_Pooled_Results
%   Prints LKB results of pooled analysis
%   Used for table 1 in paper

% Need all data
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
lkb_save_loc='Z:/elw/MATLAB/meta_analy/meta_data/lkb_parameters.mat';

pld_data = {'CGmsk','CGnki','CGrtog','CGum','CGcomb'};
protocols = {'MSK','NKI','RTOG','UMich','Comb'};

% TMP
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
pld_data = {'CGcomb'};
protocols = {'Comb'};

if isunix %on mac
    fn = strrep(fn,'Z:/elw/','/Users/elw/Documents/');
    lkb_save_loc = strrep(lkb_save_loc,'Z:/elw/','/Users/elw/Documents/');
end

%% Saved output to be read by RP_Meta_Analysis.m
% stores individual institutes lkb results, and 95% CIs
% (a,n,td50,m) X (msk, nki, rtog, um, comb)
lkb_vals = zeros(length(pld_data),4);
lkb_lcl = zeros(length(pld_data),4);
lkb_ucl = zeros(length(pld_data),4);



for i=1:length(pld_data)
    CGdata=load(fn,pld_data{i});
    CGdata = struct2cell(CGdata);
    CGdata = CGdata{1};
    
    cur_data = protocols{i};
    

    
    % find current best parameters
    llhd_mx = max(CGdata.mLymanGrid.loglikelihood(:));
    
    
    llhd_n = max(max(CGdata.mLymanGrid.loglikelihood,[],2),[],1);
    llhd_a = flipdim(squeeze(llhd_n),1);
    llhd_td50 = max(max(CGdata.mLymanGrid.loglikelihood,[],3),[],2);
    llhd_m = max(max(CGdata.mLymanGrid.loglikelihood,[],3),[],1);
    
    

     low68 = llhd_mx -0.5*(chi2inv(0.68,1));
     low95 = llhd_mx -0.5*(chi2inv(0.95,1));

    % 68% point of td50
    [td50, td50_68ci] = ConfidenceInterval(CGdata.mLymanGrid.TD50,llhd_td50, low68);
    [~, td50_95ci] = ConfidenceInterval(CGdata.mLymanGrid.TD50,llhd_td50, low95);
    td50_68ci(td50_68ci>100) = 100;
    td50_95ci(td50_95ci>100) = 100;
     
    [m, m_68ci] = ConfidenceInterval(CGdata.mLymanGrid.m,llhd_m, low68);
    [~, m_95ci] = ConfidenceInterval(CGdata.mLymanGrid.m,llhd_m, low95);
    
     [n, n_68ci] = ConfidenceInterval(CGdata.mLymanN,llhd_n, low68);
     [~, n_95ci] = ConfidenceInterval(CGdata.mLymanN,llhd_n, low95);
    % set limits
     n_68ci(n_68ci>10) = 10;
     n_95ci(n_95ci>10) = 10;
     
     [a, a_68ci] = ConfidenceInterval(CGdata.mLymanN,llhd_a, low68);    
     [~, a_95ci] = ConfidenceInterval(CGdata.mLymanN,llhd_a, low95);
     a_68ci(a_68ci<0) = 0.1;
     a_95ci(a_95ci<0) = 0.1;
    
    
    disp(['== ',cur_data,' ==']);
    disp(['n (68%) [95%]']);
    disp([num2str(n),' (',num2str(n_68ci(1),4),' - ',num2str(n_68ci(2),4),')',...
                        ' [',num2str(n_95ci(1),4),' - ',num2str(n_95ci(2),4),']']);
    disp('-');
    disp(['a (68%) [95%]']);
    disp([num2str(a),' (',num2str(a_68ci(1),4),' - ',num2str(a_68ci(2),4),')',...
                        ' [',num2str(a_95ci(1),4),' - ',num2str(a_95ci(2),4),']']);
    disp('-');
    
    disp(['TD50 (68%) [95%]']);
    disp([num2str(td50),' (',num2str(td50_68ci(1),4),' - ',num2str(td50_68ci(2),4),')',...
                        ' [',num2str(td50_95ci(1),4),' - ',num2str(td50_95ci(2),4),']']);
    disp('-');
    disp(['m (68%) [95%]']);
    disp([num2str(m),' (',num2str(m_68ci(1),4),' - ',num2str(m_68ci(2),4),')',...
                        ' [',num2str(m_95ci(1),4),' - ',num2str(m_95ci(2),4),']']);
    disp('-');

%         a, n, td50, m
       lkb_vals(i,1:4) = [a n td50 m];
       lkb_lcl(i,1:4) = [a_95ci(1) n_95ci(1) td50_95ci(1) m_95ci(1)];
       lkb_ucl(i,1:4) = [a_95ci(2) n_95ci(2) td50_95ci(2) m_95ci(2)];

end

    
    % Save LKB values for use in RP_Meta_Analysis.m, etc
    disp(['Saving LKB parameters to ',lkb_save_loc,'...']);
  save(lkb_save_loc,'lkb_vals','lkb_lcl','lkb_ucl');
end

