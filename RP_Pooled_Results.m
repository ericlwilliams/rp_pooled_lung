function [ output_args ] = RP_Pooled_Results()
%RP_Pooled_Results
%   Prints LKB results of pooled analysis
%   Used for table 1 in paper

% Need all data
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';

if isunix %on mac
    fn = strrep(fn,'Z:/elw/','/Users/elw/Documents/');
end

pld_data = {'CGnki','CGmsk','CGrtog','CGum','CGcomb'};
protocols = {'NKI','MSK','RTOG','UM','Comb'};

for i=1:length(pld_data)
    CGdata=load(fn,pld_data{i});
    CGdata = struct2cell(CGdata);
    CGdata = CGdata{1};
    
    cur_data = protocols{i};
    

    
    % find current best parameters
    llhd_mx = max(CGdata.mLymanGrid.loglikelihood(:));
    
    llhd_td50 = max(max(CGdata.mLymanGrid.loglikelihood,[],3),[],2);
    llhd_m = max(max(CGdata.mLymanGrid.loglikelihood,[],3),[],1);
    llhd_n = max(max(CGdata.mLymanGrid.loglikelihood,[],2),[],1);
    
    low68 = llhd_mx -0.5*(chi2inv(0.68,3));
    low95 = llhd_mx -0.5*(chi2inv(0.95,3));
%   low99 = mx_llhd -0.5*(chi2inv(0.99,3));

    % 68% point of td50
    [td50, td50_68ci] = ConfidenceInterval(CGdata.mLymanGrid.TD50,llhd_td50, low68);
    [~, td50_95ci] = ConfidenceInterval(CGdata.mLymanGrid.TD50,llhd_td50, low95);
    disp(['== ',cur_data,' ==']);
    disp(['TD50 (68%) [95%]']);
    disp([num2str(td50),' (',num2str(td50_68ci(1),4),' - ',num2str(td50_68ci(2),4),')',...
                        ' [',num2str(td50_95ci(1),4),' - ',num2str(td50_95ci(2),4),']']);
    
    
%     disp([num2str(td50),' - ',num2str(td50_ci)]);


    
    
%     b0sd(k) = min(CI(2)-b0(k), b0(k)-CI(1));
%     % 68% point of b1
%     [b1(k), CI] = ConfidenceInterval(OCobjs(k).mLymanGrid.m,likely2, CI68);
%     b1sd(k) = min(CI(2)-b1(k), b1(k)-CI(1));
%     % 68% point of a
%     [a(k), CI] = ConfidenceInterval(OCobjs(k).mLymanN,likely3, CI68);
%     asd(k) = min(CI(2)-a(k), a(k)-CI(1));

end

end

