function EudAnalysis_ROC_gEUD
%% Determine best gEUD a parameter using ROC analysis only
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
do_print = true;
fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';
  
fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
load(fn,'CGmsk','CGnki','CGrtog','CGum','CGcomb');
logas = -log10(CGcomb.mLymanN);
CGs = {CGmsk,CGnki,CGrtog,CGum,CGcomb};

aucs = inf(length(CGs),length(logas));
max_aucs = inf(length(CGs),1);
max_log10a = inf(length(CGs),1);

for i=1:length(CGs)

    CG = CGs{i};
   euds = [CG.mGrp.mEUD];
   
   flgcomp = ~[CG.mGrp.mFlgCensor]';
   labels=cell(size(flgcomp));
   labels(flgcomp)={'comp'};
   labels(~flgcomp)={'cens'};
   
   for j=1:length(logas)
      cur_euds = euds(j,:); 
      [~,~,~,cur_auc] = perfcurve(labels,cur_euds,'comp');
      aucs(i,j)=cur_auc;
   end
   
    %find max
    [max_aucs(i),max_auc_ind] = max(aucs(i,:));
    max_log10a(i) = logas(max_auc_ind);
end
    cur_fig=figure(1); clf reset; 
    set(cur_fig,'Position',ss_four2three);

    h_aucs=plot(logas,aucs,'LineWidth',2);hold on;
    plot(max_log10a,max_aucs,'k*','MarkerSize',14);
    
    h_aucs_lgnd = legend(h_aucs,'MSK','NKI','RTOG','UMich','Comb',...
        'Location','Best');
    set(h_aucs_lgnd,'FontSize',18);
    set(gca,'FontSize',18);
    xlabel('log_{10}(a)','FontSize',20);
    ylabel('AUC','FontSize',20);
    
    disp(['$$$']);
    disp(['Best log(a)s: ',num2str(max_log10a')])    
    disp(['$$$']);
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'h_geud_aucs'],'-png');
        disp(['Saving ',fig_loc,'h_geud_aucs.pdf...']);
    end
    
end