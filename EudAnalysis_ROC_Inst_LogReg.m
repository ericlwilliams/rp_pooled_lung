function EudAnalysis_ROC_Inst_LogReg
tic; %close all;
do_boxplots=false;
do_rtog=true;

screen_size=get(0,'ScreenSize');
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';
fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta.mat';
%fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_meta.mat';
load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');

CG_um_rtog = CGum;
CG_um_rtog = CG_um_rtog.fAddPatient(CGrtog.mGrp);

CG_um_msk = CGum;
CG_um_msk = CG_um_msk.fAddPatient(CGmsk.mGrp);

CG_um_nki = CGum;
CG_um_nki = CG_um_nki.fAddPatient(CGnki.mGrp);

CG_um_rtog_msk = CG_um_rtog;
CG_um_rtog_msk = CG_um_rtog_msk.fAddPatient(CGmsk.mGrp);

CG_um_rtog_nki = CG_um_rtog;
CG_um_rtog_nki = CG_um_rtog_nki.fAddPatient(CGnki.mGrp);

CG_um_rtog_msk_nki = CG_um_rtog_msk;
CG_um_rtog_msk_nki = CG_um_rtog_msk_nki.fAddPatient(CGnki.mGrp);


% Build models    
CGs = {CGmsk CGnki CGrtog CGum CG_um_rtog_msk_nki};
%CGs = {CG_um_rtog_msk CG_um_rtog_nki CG_um_rtog_msk_nki};
%groups = {'UMich+RTOG+MSK' 'UMich+RTOG+NKI' 'UMich+RTOG+MSK+NKI'};
groups = {'MSK' 'NKI' 'RTOG' 'UMich' 'MSK+NKI+RTOG+UMich'};

colors = {'r' 'b' 'g' 'm' 'k'};
cur_fig=figure;clf reset;
aucs = zeros(length(CGs),1);
pvals = zeros(length(CGs),1);
set(gcf,'Position',ss_four2three);
plot([0 1],[0 1],'--');hold on;        
for i=1:length(CGs)
    CG = CGs{i};
    st = [CG.mLogisticRegressionMat];
    dpf = [st.dev]; % deviations
    [~,loc] = min(dpf); % the minimum deviation corresponds to the maximum likelihood
    disp(['best log10(a) of Logistic Regression of exact gEUD is: ',num2str(-CG.mLymanN(loc))]);
    loga = -CG.mLymanN(loc);
    
    [~,loc] = min(abs(CG.mLymanN+loga)); % the n whose corresponding responding function will be ploted
    disp(['the log10(a) in responding curve is: ',num2str(loga)]);

    st = CG.mLogisticRegressionMat(loc); % the fitting result of that n
    euds = [CG.mGrp.mEUD]; euds = euds(loc,:); % the gEUDs of that n

    pvalue = st.stats.p;
    pval = pvalue(2); % the p-value corresponding to gEUD

    %doses = (0:max(euds))'; % doses (gEUD) whose RP probability will be computed
    %dose = [CG.mGrp.mDoseBins_LQ]';
        
    
    [rpb,~,~] = glmval(st.b, euds,'logit',st.stats); % the responding function values at doses
    disp(['the beta values are: ',num2str(st.b')]);
    
    flgcomp = ~[CG.mGrp.mFlgCensor]';
    labels=cell(size(flgcomp));
    labels(flgcomp)={'comp'};
    labels(~flgcomp)={'cens'};
    %roc(rpb,flgcomp);     
    %[tpr,fpr,thresholds]=roc(flgcomp,rpb);
    % plotroc(flgcomp,rpb,'test');
    [X,Y,~,AUC] = perfcurve(labels,rpb,'comp');
    aucs(i)=AUC;
    
    %standard error of area
    Area=AUC;
    lu = sum(flgcomp);
    lh = sum(~flgcomp);
    Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
    V=(Area*(1-Area)+(lu-1)*(Q1-Area2)+(lh-1)*(Q2-Area2))/(lu*lh);
    Serror=realsqrt(V);
    SAUC=(Area-0.5)/Serror; %standardized area
    pvals(i)=1-0.5*erfc(-SAUC/realsqrt(2)); %p-value

    h(i)=plot(X,Y,colors{i},'LineWidth',2);
    xlabel('False positive rate (1-Specificity)','FontSize',14);
    ylabel('True positive rate (Sensitivity)','FontSize',14)
    title('ROC for classification by logistic regression','FontSize',16)
end
    set(gca,'FontSize',12);
    h_lgnd=legend(h,...
       [groups{1},10,'AUC: ',num2str(aucs(1),3)],...
       [groups{2},10,'AUC: ',num2str(aucs(2),3)],...
       [groups{3},10,'AUC: ',num2str(aucs(3),3)],...
       [groups{4},10,'AUC: ',num2str(aucs(4),3)],...
       [groups{5},10,'AUC: ',num2str(aucs(5),3)],...
        'Location','Best');
    set(h_lgnd,'FontSize',14);
    
       
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'logreg_roc'],'-pdf');
                disp(['Saving ',fig_loc,'logreg_roc.pdf...']);
     


end