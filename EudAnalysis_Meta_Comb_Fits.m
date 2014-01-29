function EudAnalysis_Meta_Comb_Fits
tic; %close all;

screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

    do_print = true;
    fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';
    % load results
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_RTOG_EUD_meta_lt6m_vnorm.mat';
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_RTOG_Homo_EUD_meta.mat';
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_RTOG_phys_a2b_Hetero_EUD_meta_lt6m.mat';
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_ntd_a2b_Hetero_EUD_meta_lt6m.mat';
    
    
    
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_meta.mat';
    %fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_fine_meta.mat';
    %fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_med_meta.mat';
   
    
    
    %fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_crs_meta.mat';
   
    %% sskipefilm1
    %fn = 'C:\Documents and Settings\williae1\meta_data\MSK_NKI_UMich_RTOG_EUD_crs_meta.mat';
    
    %fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_med_EUD_fine_meta_comb.mat';
    fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';

    %     do_a2b=true;
%     a2b='phys';
%     

    load(fn,'CGcomb');
    
% prepare
    LymanN = log10(CGcomb.mLymanN);

    CGcomb.mLymanN = LymanN;
    
    cm2 = colormap(jet(10));

% lyman analysis p-values
if 1
    ctr=1;
    for i=length(LymanN):-1:1,
        cur_loga=-LymanN(i);
        cur_fig=figure(i); clf reset; colormap(cm2);
        set(gcf,'Position',2.*ss_four2three);

%         subplot(2,4,[1:2]);
%         h_llhd_comb=CGcomb.fLymanGridExactFig_a_loglikelihood('ks--',2);
%         line([cur_loga cur_loga],ylim,'LineWidth',3,'Color','r');
%         
%         subplot(2,4,[3:4]);
%         h_pval_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
%         line([cur_loga cur_loga],ylim,'LineWidth',3,'Color','r');
%     
%   
%         subplot(2,4,[6:7])
%         [loga,pval] = CGcomb.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'r',1);
%         CGcomb.fComplicationObservedFig_EUD(loga,4,'r*',2);
%         %xlim([0 90]);
%         loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
%         str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
%         text(0.65,0.85,str,'FontSize',24,'Units','normalized');
%         %text(2,0.6,str,'FontSize',18);
%         title('MSK + NKI + RTOG + UMich','FontSize',24);    

       subplot(2,4,[2:3]);
       h_llhd_comb=CGcomb.fLymanGridExactFig_a_loglikelihood('ks--',2);
         line([cur_loga cur_loga],ylim,'LineWidth',3,'Color','r');
                 
         subplot(2,4,[6:7])
        [loga,pval] = CGcomb.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'r',1);
        CGcomb.fComplicationObservedFig_EUD(loga,4,'r*',2);
        %xlim([0 90]);
        loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
        str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
        text(0.65,0.85,str,'FontSize',24,'Units','normalized');
        %text(2,0.6,str,'FontSize',18);
        title('MSK + NKI + RTOG + UMich','FontSize',24); 
        
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'zoom_rp_lkb_comb2-',num2str(ctr)],'-pdf');
        disp(['Saving ',fig_loc,'zoom_rp_lkb_comb2-',num2str(ctr),'.pdf...']);
    end
    ctr=ctr+1;
    end
end    





%# Logistic Regression Response function for all loga
% 4x4 with subplots
% if 1
%         for i=1:length(LymanN),
%             cur_loga=-LymanN(i);
%             cur_fig=figure(i);clf reset;colormap(cm2);
%             set(cur_fig,'Position',[0 0 screen_size(3) screen_size(4)]);
%             hold off;
%             subplot(2,2,2)
%             [loga,pval] = CGrtog.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'g',1);
%             CGrtog.fComplicationObservedFig_EUD(loga,4,'g*',1);
%             if loga~=cur_loga,
%                 disp(['loga ~= cur_loga!']);
%             end
%             loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
%             str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
%             text(0.75,0.85,str,'FontSize',16,'Units','normalized');
%             %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
%             %text(2,0.6,str,'FontSize',18);
%             title('RTOG','FontSize',24);
%             
%             subplot(2,2,1)
%             [loga,pval] = CGmsk.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'r',1);
%             CGmsk.fComplicationObservedFig_EUD(loga,4,'r*',1);
%             loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
%             str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
%             text(0.75,0.85,str,'FontSize',16,'Units','normalized');
%             %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
%             %text(2,0.6,str,'FontSize',18);
%             title('MSK','FontSize',24);
%             
%             subplot(2,2,4)
%             [loga,pval] = CGnki.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'b',1);
%             CGnki.fComplicationObservedFig_EUD(loga,4,'b*',1);
%             loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
%             str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));            
%             text(0.75,0.85,str,'FontSize',16,'Units','normalized');
%             %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
%             %text(2,0.6,str,'FontSize',18);
%             title('NKI','FontSize',24);
%             
%             subplot(2,2,3)
%             [loga,pval] = CGum.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'m',1);
%             CGum.fComplicationObservedFig_EUD(loga,4,'m*',1);
%             loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
%             str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
%             text(0.75,0.85,str,'FontSize',16,'Units','normalized');
%             %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
%             %text(2,0.6,str,'FontSize',18);
%             title('UMich','FontSize',24);
%             
%             if do_print,
%                 set(cur_fig,'Color','w');
%                 export_fig(cur_fig,[fig_loc,'logreg_responses_frange',...
%                     '-',num2str(i)],'-pdf');
%                 disp(['Saving ',fig_loc,'logreg_responses_frange',...
%                     '-',num2str(i),'.pdf...']);
%             end;
%         end;
% end;





toc;