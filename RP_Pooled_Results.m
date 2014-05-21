function  RP_Pooled_Results()
%RP_Pooled_Results
%   Prints LKB results of pooled analysis
%   Used for table 1 in paper

% Need all data
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
pld_data = {'CGnki','CGmsk','CGrtog','CGum','CGcomb'};

if isunix %on mac
    fn = strrep(fn,'Z:/elw/','/Users/elw/Documents/');
end


protocols = {'NKI','MSK','RTOG','UMich','Comb'};

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
    
    
    % Fan's method (assumes df=1), underestimate confidence intervals     
     low68 = llhd_mx-0.5* (1);
     low95 = llhd_mx-0.5* (1.96^2);

%      low68 = llhd_mx -0.5*(chi2inv(0.68,3));
%      low95 = llhd_mx -0.5*(chi2inv(0.95,3));

    % 68% point of td50
    [td50, td50_68ci] = ConfidenceInterval(CGdata.mLymanGrid.TD50,llhd_td50, low68);
    [~, td50_95ci] = ConfidenceInterval(CGdata.mLymanGrid.TD50,llhd_td50, low95);
    
    [m, m_68ci] = ConfidenceInterval(CGdata.mLymanGrid.m,llhd_m, low68);
    [~, m_95ci] = ConfidenceInterval(CGdata.mLymanGrid.m,llhd_m, low95);
    
     [n, n_68ci] = ConfidenceInterval(CGdata.mLymanN,llhd_n, low68);
     [~, n_95ci] = ConfidenceInterval(CGdata.mLymanN,llhd_n, low95);

     [a, a_68ci] = ConfidenceInterval(CGdata.mLymanN,llhd_a, low68);    
     [~, a_95ci] = ConfidenceInterval(CGdata.mLymanN,llhd_a, low95);
        
    
    
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

            

    
end

end

