function RP_Scratch

% Only loading combined likelihoods
fname='Z:/elw/MATLAB/meta_analy/meta_data/lkb_comb_llhds.mat';

    
if isunix %on mac
    fname = strrep(fname,'Z:/elw/','/Users/elw/Documents/');
end
load(fname,'lkb_llhds','lkb_n','lkb_td50','lkb_m');


% use mx_llhd_llhd to find contours when slicing at meta points
[llhd_mx, mx_llhd_ind] = max(lkb_llhds(:));
[comb_td50_idx, comb_m_idx, comb_n_idx] = ind2sub(size(lkb_llhds),mx_llhd_ind);

% Define confidence intervals
% Fan's method (assumes df=1), underestimate confidence intervals


 %% Method 1 <- 'full' cis?
 
 llhd_n = max(max(lkb_llhds,[],2),[],1);
 llhd_a = flipdim(squeeze(llhd_n),1);
 llhd_td50 = max(max(lkb_llhds,[],3),[],2);
 llhd_m = max(max(lkb_llhds,[],3),[],1);

 low68 = llhd_mx-0.5* (1);
 low95 = llhd_mx-0.5* (1.96^2); 
 
 [td50, td50_68ci] = ConfidenceInterval(lkb_td50,llhd_td50, low68);
 [~, td50_95ci] = ConfidenceInterval(lkb_td50,llhd_td50, low95);
 td50_68ci(td50_68ci>100) = 100;
 td50_95ci(td50_95ci>100) = 100;
 
 [m, m_68ci] = ConfidenceInterval(lkb_m,llhd_m, low68);
 [~, m_95ci] = ConfidenceInterval(lkb_m,llhd_m, low95);
 
 [n, n_68ci] = ConfidenceInterval(lkb_n,llhd_n, low68);
 [~, n_95ci] = ConfidenceInterval(lkb_n,llhd_n, low95);
 % set limits
 n_68ci(n_68ci>10) = 10;
 n_95ci(n_95ci>10) = 10;
 
 [a, a_68ci] = ConfidenceInterval(lkb_n,llhd_a, low68);
 [~, a_95ci] = ConfidenceInterval(lkb_n,llhd_a, low95);
 a_68ci(a_68ci<0) = 0.1;
 a_95ci(a_95ci<0) = 0.1;

 
  disp(['== Method 1 ==']);
    disp(['n (68%) [95%] <95% diff>']);
    disp([num2str(n),' (',num2str(n_68ci(1),4),' - ',num2str(n_68ci(2),4),')',...
                        ' [',num2str(n_95ci(1),4),' - ',num2str(n_95ci(2),4),']',...
                        ' <',num2str(abs(n_95ci(2)-n_95ci(1)),4),'>']);
    disp('-');
    disp(['a (68%) [95%]']);
    disp([num2str(a),' (',num2str(a_68ci(1),4),' - ',num2str(a_68ci(2),4),')',...
                        ' [',num2str(a_95ci(1),4),' - ',num2str(a_95ci(2),4),']',...
                        ' <',num2str(abs(a_95ci(2)-a_95ci(1)),4),'>']);
    disp('-');
    
    disp(['TD50 (68%) [95%]']);
    disp([num2str(td50),' (',num2str(td50_68ci(1),4),' - ',num2str(td50_68ci(2),4),')',...
                        ' [',num2str(td50_95ci(1),4),' - ',num2str(td50_95ci(2),4),']',...
                        ' <',num2str(abs(td50_95ci(2)-td50_95ci(1)),4),'>']);
    disp('-');
    disp(['m (68%) [95%]']);
    disp([num2str(m),' (',num2str(m_68ci(1),4),' - ',num2str(m_68ci(2),4),')',...
                        ' [',num2str(m_95ci(1),4),' - ',num2str(m_95ci(2),4),']',...
                        ' <',num2str(abs(m_95ci(2)-m_95ci(1)),4),'>']);
    disp('-');
 
 %% Method 2
 % split at best n, but use 2df technique from there
 td50_m_llhds = lkb_llhds(:,:,comb_n_idx);

 df=1;
% low68 = llhd_mx -0.5*(chi2inv(0.68,df));
% low95 = llhd_mx -0.5*(chi2inv(0.95,df)); 
  
ctr_matrix=contourc(lkb_td50,lkb_m,td50_m_llhds',[low95 low68]);
end_95ci_idx = ctr_matrix(2,1);
ctr_mtx_95ci = ctr_matrix(:,2:end_95ci_idx+1);
ctr_mtx_68ci = ctr_matrix(:,end_95ci_idx+3:end);

ctr_td50_68ci = [min(ctr_mtx_68ci(1,:)) max(ctr_mtx_68ci(1,:))];
ctr_td50_95ci = [min(ctr_mtx_95ci(1,:)) max(ctr_mtx_95ci(1,:))];

ctr_m_68ci = [min(ctr_mtx_68ci(2,:)) max(ctr_mtx_68ci(2,:))];
ctr_m_95ci = [min(ctr_mtx_95ci(2,:)) max(ctr_mtx_95ci(2,:))];


disp(['== Method 2 (df = ',num2str(df),') ==']);

 
    disp(['TD50 (68%) [95%]']);
    disp([num2str(lkb_td50(comb_td50_idx)),' (',num2str(ctr_td50_68ci(1),4),' - ',num2str(ctr_td50_68ci(2),4),')',...
                        ' [',num2str(ctr_td50_95ci(1),4),' - ',num2str(ctr_td50_95ci(2),4),']',...
                        ' <',num2str(abs(ctr_td50_95ci(2)-ctr_td50_95ci(1)),4),'>']);
    disp('-');
    disp(['m (68%) [95%]']);
    disp([num2str(lkb_m(comb_m_idx)),' (',num2str(ctr_m_68ci(1),4),' - ',num2str(ctr_m_68ci(2),4),')',...
                        ' [',num2str(ctr_m_95ci(1),4),' - ',num2str(ctr_m_95ci(2),4),']',...
                        ' <',num2str(abs(ctr_m_95ci(2)-ctr_m_95ci(1)),4),'>']);
    disp('-');

    
%% Method 3, cut through slice, then maximize like Method 1
    
m3_llhd_m = max(td50_m_llhds,[],1);
m3_llhd_td50 = max(td50_m_llhds,[],2);

low68 = llhd_mx-0.5* (1);
 low95 = llhd_mx-0.5* (1.96^2);

 [td50, td50_68ci] = ConfidenceInterval(lkb_td50,m3_llhd_td50, low68);
 [~, td50_95ci] = ConfidenceInterval(lkb_td50,m3_llhd_td50, low95);
 td50_68ci(td50_68ci>100) = 100;
 td50_95ci(td50_95ci>100) = 100;
 
 [m, m_68ci] = ConfidenceInterval(lkb_m,m3_llhd_m, low68);
 [~, m_95ci] = ConfidenceInterval(lkb_m,m3_llhd_m, low95);
 
 
  disp(['== Method 3 ==']);
    disp(['TD50 (68%) [95%]']);
    disp([num2str(td50),' (',num2str(td50_68ci(1),4),' - ',num2str(td50_68ci(2),4),')',...
                        ' [',num2str(td50_95ci(1),4),' - ',num2str(td50_95ci(2),4),']',...
                        ' <',num2str(abs(td50_95ci(2)-td50_95ci(1)),4),'>']);
    disp('-');
    disp(['m (68%) [95%]']);
    disp([num2str(m),' (',num2str(m_68ci(1),4),' - ',num2str(m_68ci(2),4),')',...
                        ' [',num2str(m_95ci(1),4),' - ',num2str(m_95ci(2),4),']',...
                        ' <',num2str(abs(m_95ci(2)-m_95ci(1)),4),'>']);
 
disp('');
end

