function EudAnalysis_Meta_disp
tic; %close all;

%% random effects meta analysis results [val,95% lo, hi]
re_td50 = [21.88,15.9,27.8];
re_log10a = [-0.0862,-0.69897,0.15533];
re_m = [0.29,0.13,0.44];

%% fixed effects 
fe_td50 = [21.88,15.9,27.8];
fe_log10a = [-0.0506,-0.2924,0.1038];
fe_m = [0.21,0.14,0.28];

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
    %fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_med_EUD_fine_meta.mat';
   
    
    
    %fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
    %fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine
    %_meta.mat';
    
    %% COMB only
    %fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
    %is_comb_only=true;
    fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
    is_comb_only=false;
    
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_crs_meta.mat';
   
    %% sskipefilm1
    %fn = 'C:\Documents and Settings\williae1\meta_data\MSK_NKI_UMich_RTOG_EUD_crs_meta.mat';
    
    %fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';
    
    %fn ='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_med_meta.mat';
    
    
    %fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';

    %     do_a2b=true;
%     a2b='phys';
%     
    if isunix
        fn = strrep(fn, 'G:', '/media/SKI_G');
        fn_msk = strrep(fn_msk, 'G:', '/media/SKI_G');
    end
    
    if is_comb_only
        load(fn,'CGcomb');
    else
        load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');
    end
    
    numreg=1;
    Regions = {'Whole'};
    
% prepare
    cm1 = colormap(jet(300)); cm1=cm1(1:256,:); %cm1(end,:) = 0.5;
    cm2 = colormap(jet(10));
    ticksz = 10;
    fontsz = 10;
     % log10(n) correction for figures
    
    if ~is_comb_only,
        LymanN = log10(CGmsk.mLymanN);
        CGmsk.mLymanN = LymanN;
        CGnki.mLymanN = LymanN;
        CGum.mLymanN = LymanN;
        CGrtog.mLymanN = LymanN;
        CGcomb.mLymanN = LymanN;
        CGtmp.mLymanN = LymanN;
    else
        LymanN = log10(CGcomb.mLymanN);
        CGcomb.mLymanN = LymanN;
    end
    
    
    %% plot log-likelihood isosurface
    if 1
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        
        CGcomb.fLymanGridExactFig_IsoSurface();
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_3d'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_3d.png...']);
        end
        
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_IsoSurface(1);
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_3d_p1'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_3d_p1.png...']);
        end
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_IsoSurface(2);
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_3d_p2'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_3d_p2.png...']);
        end
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_IsoSurface(3);
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_3d_p3'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_3d_p3.png...']);
         end
        
    end
    
% plot DVHs -> EudAnalysis_DVHs.m
if 0
   CGnki.fDVHCurvesSummary_DVH();
%    
%     disp('RTOG DVH curves:');
%   
%     figure(1); clf reset; hold on; % grid on;
%     f = [CGnki.mGrp.mFlgCensor];
%     % DVHs of censored patients
%     g = find(f);
%     for k = 1:length(g)
%         %plot(CGnki.mGrp(g(k)).mDoseBins_org, CGnki.mGrp(g(k)).mVolCum);
%         plot(CGnki.mGrp(g(k)).mDoseBins_org, CGnki.mGrp(g(k)).mVolDiff);
%     end
%     % DVHs of complicated patients
%     g = find(~f);
%     for k = 1:length(g)
%         %plot(CGnki.mGrp(g(k)).mDoseBins_org, CGnki.mGrp(g(k)).mVolCum,'r');
%         plot(CGnki.mGrp(g(k)).mDoseBins_org, CGnki.mGrp(g(k)).mVolDiff,'r');
%     end
%     ylim([0,0.5]);
%     set(gca,'xminortick','on','yminortick','on');
% %     
%     figure(1); clf reset; hold on; % grid on;
%     f = [CGrtog.mGrp.mFlgCensor];
%     % DVHs of censored patients
%     for k=1:length(f),
%         if f(k), %# cens
%             dvh_color = 'b';
%         else
%             dvh_color = 'r';
%         end
%         dosebins = CGrtog.mGrp(k).mDoseBins_org;
%         volcum = CGrtog.mGrp(k).mVolCum;
%         plot(dosebins(1:end-1), volcum(1:end-1),dvh_color);
%     end
%     set(gca,'xminortick','on','yminortick','on');
%     xlabel('Dose [Gy]','FontSize',15);
%     ylabel('Volume Fraction','FontSize',15);
%     
%     figure(2); clf reset; hold on; % grid on;
%     f = [CGrtog.mGrp.mFlgCensor];
%     % DVHs of censored patients
%     for k=1:length(f),
%         if f(k), %# cens
%             dvh_color = 'b';
%         else
%             dvh_color = 'r';
%         end
%         dosebins = CGrtog.mGrp(k).mDoseBins_org;
%         voldiff = CGrtog.mGrp(k).mVolDiff;
%         plot(dosebins(1:end-1), voldiff(1:end-1),dvh_color);
%     end
%     set(gca,'xminortick','on','yminortick','on');
%     xlabel('Dose [Gy]','FontSize',15);
%     ylabel('Volume Fraction','FontSize',15);
%     
end

% atlases
if 0
    
    if ~is_comb_only,
        %disp(' '); disp(Regions{k});
        figure(1); clf reset;
        CGmsk.fAtlasCompactFig_EUD(fontsz,ticksz);
        xlabel(''); ylabel('');
        figure(2); clf reset;
        CGnki.fAtlasCompactFig_EUD(fontsz,ticksz);
        xlabel(''); ylabel('');
        figure(3); clf reset;
        CGrtog.fAtlasCompactFig_EUD(fontsz,ticksz);
        xlabel(''); ylabel('');
    end
    figure(4); clf reset;
    CGcomb.fAtlasCompactFig_EUD(fontsz,ticksz);
    xlabel(''); ylabel('');
end

% dvh atlases - moved to scripts/GenerateDvhAtlas.m
if 0
    
    if ~is_comb_only,
        %disp(' '); disp(Regions{k});
        figure(1); clf reset;
        CGmsk.fCrudeAtlas_DVH(-1);
        xlabel(''); ylabel('');
        figure(2); clf reset;
        CGnki.fCrudeAtlas_DVH(-1);
        xlabel(''); ylabel('');
        figure(3); clf reset;
        CGrtog.fCrudeAtlas_DVH(-1);
        xlabel(''); ylabel('');
    end
    figure(4); clf reset;
    CGcomb.fCrudeAtlas_DVH(-1); %nfx = -1, all fx schemes
    xlabel(''); ylabel('');
end

% probability of having >=20% RP rate
if 0
  
  %disp(' '); disp(Regions{k});
  cur_fig=figure(1); clf reset;
  %set(gcf,'Position',ss_four2three);
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGmsk.fProbabilityFig_EUD('');
  %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
  if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'msk_rp20pct'],'-pdf');
   disp(['Saving ',fig_loc,'msk_rp20pct.pdf...']);
  end
  xlabel(''); ylabel('');

  cur_fig=figure(2); clf reset;
  %set(gcf,'Position',ss_four2three);
  
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGnki.fProbabilityFig_EUD('');
  %if do_print, print_fig(gcf,fig_loc,'nki_rp20pct','pdf');end;
  if do_print,
      set(cur_fig,'Color','w');
      export_fig(cur_fig,[fig_loc,'nki_rp20pct'],'-pdf');
      disp(['Saving ',fig_loc,'nki_rp20pct.pdf...']);
  end
  xlabel(''); ylabel('');
  
  cur_fig=figure(3); clf reset;
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  %set(gcf,'Position',ss_four2three);
   CGum.fProbabilityFig_EUD('');
  %if do_print, print_fig(gcf,fig_loc,'um_rp20pct','pdf');end;
  if do_print,
      set(cur_fig,'Color','w');
      export_fig(cur_fig,[fig_loc,'um_rp20pct'],'-pdf');
      disp(['Saving ',fig_loc,'um_rp20pct.pdf...']);
  end
  xlabel(''); ylabel('');
  
  cur_fig=figure(4); clf reset;
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  %set(gcf,'Position',ss_four2three);
  CGrtog.fProbabilityFig_EUD('');
  %if do_print, print_fig(gcf,fig_loc,'rtog_rp20pct','pdf');end;
  if do_print,
      set(cur_fig,'Color','w');
      export_fig(cur_fig,[fig_loc,'rtog_rp20pct'],'-pdf');
      disp(['Saving ',fig_loc,'rtog_rp20pct.pdf...']);
  end
  xlabel(''); ylabel('');
  
  cur_fig=figure(5); clf reset;
  set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  %set(gcf,'Position',ss_four2three);
  CGcomb.fProbabilityFig_EUD('');
  %if do_print, print_fig(gcf,fig_loc,'comb_rp20pct','png');end;
  if do_print,
      set(cur_fig,'Color','w');
      export_fig(cur_fig,[fig_loc,'comb_rp20pct'],'-pdf');
      disp(['Saving ',fig_loc,'comb_rp20pct.pdf...']);
  end
  xlabel(''); ylabel('');

  cur_fig=figure(6);clf reset;
  %set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
  CGtmp = {CGmsk CGrtog CGum CGnki};
  CGtitles = {'MSK' 'RTOG' 'UMich' 'NKI'};
  for k=1:4,
    if k==1, set(gcf,'Position',2.*ss_four2three);end;
    subplot(2,2,k);
    hold on;
    %title(CGtitles{k},'FontSize',16);
    CGtmp{k}.fProbabilityFig_EUD(CGtitles{k});
    
  end
    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'rp20pct'],'-pdf');
        disp(['Saving ',fig_loc,'rp20pct.pdf...']);
    end;
  
  
end
% low 68% confidence
if 0
   
    cur_fig=figure(1); clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    CGmsk.fLow68pctConfidenceFig_EUD();
    %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'msk_68cl'],'-pdf');
        disp(['Saving ',fig_loc,'msk_68cl.pdf...']);
    end
    xlabel(''); ylabel('');
    
    cur_fig=figure(2); clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    CGnki.fLow68pctConfidenceFig_EUD();
    %if do_print, print_fig(gcf,fig_loc,'msk_rp20pct','pdf');end;
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'nki_68cl'],'-pdf');
        disp(['Saving ',fig_loc,'nki_68cl.pdf...']);
    end
    xlabel(''); ylabel('');
    
    
     cur_fig=figure(3); clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    CGum.fLow68pctConfidenceFig_EUD();
    %if do_print, print_fig(gcf,fig_loc,'msk_68cl','pdf');end;
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'um_68cl'],'-pdf');
        disp(['Saving ',fig_loc,'um_68cl.pdf...']);
    end
    xlabel(''); ylabel('');
    
      cur_fig=figure(4); clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    CGrtog.fLow68pctConfidenceFig_EUD();
    %if do_print, print_fig(gcf,fig_loc,'msk_68cl','pdf');end;
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'rtog_68cl'],'-pdf');
        disp(['Saving ',fig_loc,'rtog_68cl.pdf...']);
    end
    xlabel(''); ylabel('');
    
    cur_fig=figure(5); clf reset;
    set(cur_fig,'Position',[0 0 screen_size(3)/2 screen_size(4)/2]);
    CGcomb.fLow68pctConfidenceFig_EUD();
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'comb_68cl'],'-pdf');
        disp(['Saving ',fig_loc,'comb_68cl.pdf...']);
    end
    xlabel(''); ylabel('');
    
    %     %disp(' '); disp(Regions{k});
    %     figure(1); clf reset;
%     set(gcf,'Position',ss_four2three);
%     CGmsk.fLow68pctConfidenceFig_EUD();
%     if do_print, print_fig(gcf,fig_loc,'msk_68cl','png');end;
%     xlabel(''); ylabel('');
%     
%     figure(2); clf reset;
%     set(gcf,'Position',ss_four2three);
%     CGnki.fLow68pctConfidenceFig_EUD();
%     if do_print, print_fig(gcf,fig_loc,'nki_68cl','png');end;
%     xlabel(''); ylabel('');
%     
%     figure(3); clf reset;
%     set(gcf,'Position',ss_four2three);
%     CGum.fLow68pctConfidenceFig_EUD();
%     if do_print, print_fig(gcf,fig_loc,'um_68cl','png');end;
%     xlabel(''); ylabel('');
%     
%     figure(4); clf reset;
%     set(gcf,'Position',ss_four2three);
%     CGrtog.fLow68pctConfidenceFig_EUD();
%     if do_print, print_fig(gcf,fig_loc,'rtog_68cl','png');end;
%     xlabel(''); ylabel('');
%     
%     figure(5); clf reset;
% set(gcf,'Position',ss_four2three);
%     CGcomb.fLow68pctConfidenceFig_EUD();
%     if do_print, print_fig(gcf,fig_loc,'comb_68cl','png');end;
%     xlabel(''); ylabel('');

end

% G-value of goodness of fit
if 0
    disp(['p-value of goodness of fit: Lyman Model Whole lung: ',10,...
            'MSK: ',num2str(CGmsk.mLymanGoodnessOfFitSim.p_value),10,...
            'NKI: ',num2str(CGnki.mLymanGoodnessOfFitSim.p_value),10,...
            'UMich: ',num2str(CGum.mLymanGoodnessOfFitSim.p_value),10,...
            'RTOG: ',num2str(CGrtog.mLymanGoodnessOfFitSim.p_value),10,...
            'COMB: ',num2str(CGcomb.mLymanGoodnessOfFitSim.p_value)]);

        disp(['p-value of goodness of fit: Logistic Model Whole lung: ',10,...
            'MSK: ',num2str(CGmsk.mLogisticRegressionGoodnessOfFitSim.p_value),10,...
            'NKI: ',num2str(CGnki.mLogisticRegressionGoodnessOfFitSim.p_value),10,...
            'UMich: ',num2str(CGum.mLogisticRegressionGoodnessOfFitSim.p_value),10,...
            'RTOG: ',num2str(CGrtog.mLogisticRegressionGoodnessOfFitSim.p_value),10,...
            'COMB: ',num2str(CGcomb.mLogisticRegressionGoodnessOfFitSim.p_value)]);


        CGcomb.fLymanGoodnessOfFitSimulationExact_plot_EUD();
        
%     disp(['p-value of goodness of fit: Logistic Regression Whole lung: '...
%             num2str([CGmsk.mLogisticRegressionGoodnessOfFitSim.p_value,...
%             CGnki.mLogisticRegressionGoodnessOfFitSim.p_value,...
%             CGum.mLogisticRegressionGoodnessOfFitSim.p_value,...
%             CGrtog.mLogisticRegressionGoodnessOfFitSim.p_value,...
%             CGcomb.mLogisticRegressionGoodnessOfFitSim.p_value])]);
    
end

% Lyman model
% maps - TD50 & m  toprint d
if 0
          
    %disp(' '); disp(Regions{k});
      if ~is_comb_only, 
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        
        CGmsk.fLymanGridExactFig_TD50_m_EUD();
        set(gca,'YLim',[0,1.6]);
        set(gca,'box','on');
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_cntr_ll_td50_m'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_cntr_ll_td50_m.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        
        CGnki.fLymanGridExactFig_TD50_m_EUD();
        set(gca,'YLim',[0,1.6]);
        set(gca,'box','on');

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_cntr_ll_td50_m'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_cntr_ll_td50_m.png...']);
        end
        
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanGridExactFig_TD50_m_EUD();
        set(gca,'YLim',[0,1.6]);
        set(gca,'box','on');

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_cntr_ll_td50_m'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_cntr_ll_td50_m.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanGridExactFig_TD50_m_EUD();
        set(gca,'YLim',[0,1.6]);
        set(gca,'box','on');
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_cntr_ll_td50_m'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_cntr_ll_td50_m.png...']);
        end
      end %end comb only
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_TD50_m_EUD();
        set(gca,'YLim',[0,1.6]);
        set(gca,'box','on');

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_cntr_ll_td50_m'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_cntr_ll_td50_m.png...']);
        end
    
        
        %% Zoom with meta analyses

        cur_fig=figure(6); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_TD50_m_EUD();hold on;
        ylim([0.1 0.6])
        xlim([15 30])
        
        plot(re_td50(1),re_m(1),'mx','MarkerSize',15,'LineWidth',2);
        h_rem=errorbar(re_td50(1),re_m(1),re_m(1)-re_m(2),re_m(3)-re_m(1),'m.-','LineWidth',2);
        h_rem_ex=errorbar_x(re_td50(1),re_m(1),re_td50(1)-re_td50(2),re_td50(3)-re_td50(1),'m.-');
        set(h_rem_ex,'LineWidth',2);
        
        plot(fe_td50(1),fe_m(1),'gx','MarkerSize',15,'LineWidth',2);
        h_fem=errorbar(fe_td50(1),fe_m(1),fe_m(1)-fe_m(2),fe_m(3)-fe_m(1),'g.-','LineWidth',2);
        h_fem_ex=errorbar_x(fe_td50(1),fe_m(1),fe_td50(1)-fe_td50(2),fe_td50(3)-fe_td50(1),'g.-');
        
        set(h_fem_ex,'LineWidth',2);
       
        set(gca,'box','on');
        h_lgnd=legend([h_fem h_rem],'Fixed Effects (95% CI)','Random Effects (95% CI)','Location','Best');
        set(h_lgnd,'FontSize',18);
        % errorbar(medianeud,prob,max(0,prob-betainv16),max(0,betainv84-prob),strMarker,'LineWidth',lw);
        %  errorbar_x(medianeud,prob,(medianeud-binlow),(binhigh-medianeud)
        %  ,strMarker);
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_td50_m_with_meta'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_td50_m_with_meta.png...']);
        end
    
      
end
    
%profile vs m - TD50 & a
if 0
          
    %disp(' '); disp(Regions{k});
    
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGmsk.fLymanGridExactFig_Profile_TD50_a_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_prof_td50_a'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_prof_td50_a.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGnki.fLymanGridExactFig_Profile_TD50_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_prof_td50_a'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_prof_td50_a.png...']);
        end
        
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanGridExactFig_Profile_TD50_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_prof_td50_a'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_prof_td50_a.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanGridExactFig_Profile_TD50_a_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_prof_td50_a'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_prof_td50_a.png...']);
        end
        
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_Profile_TD50_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_prof_td50_a'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_prof_td50_a.png...']);
        end
    
end

%profile vs TD50 - m & a
if 0
          
    %disp(' '); disp(Regions{k});
       
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGmsk.fLymanGridExactFig_Profile_m_a_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_prof_m_a'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_prof_m_a.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGnki.fLymanGridExactFig_Profile_m_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_prof_m_a'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_prof_m_a.png...']);
        end
        
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanGridExactFig_Profile_m_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_prof_m_a'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_prof_m_a.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanGridExactFig_Profile_m_a_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_prof_m_a'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_prof_m_a.png...']);
        end
        
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_Profile_m_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_prof_m_a'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_prof_m_a.png...']);
        end
    
end


%profile vs a - td50 & m
if 0
          
    %disp(' '); disp(Regions{k});
       
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGmsk.fLymanGridExactFig_Profile_td50_m_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_prof_td50_m'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_prof_td50_m.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGnki.fLymanGridExactFig_Profile_td50_m_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_prof_td50_m'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_prof_td50_m.png...']);
        end
        
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanGridExactFig_Profile_td50_m_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_prof_td50_m'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_prof_td50_m.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanGridExactFig_Profile_td50_m_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_prof_td50_m'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_prof_td50_m.png...']);
        end
        
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_Profile_td50_m_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_prof_td50_m'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_prof_td50_m.png...']);
        end
    
end



% response vs m
if 0
          
    %disp(' '); disp(Regions{k});
       
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        loga=CGmsk.fLymanResponseFig_m_EUD(0.4);
        CGmsk.fComplicationObservedFig_EUD(loga,4,'r*',1);
         
         
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_resp_m'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_resp_m.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        loga=CGnki.fLymanResponseFig_m_EUD(0.36);
        CGnki.fComplicationObservedFig_EUD(loga,4,'r*',1);
        xlim([0 60]);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_resp_m036'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_resp_m036.png...']);
        end
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        loga=CGnki.fLymanResponseFig_m_EUD(0.38);
        CGnki.fComplicationObservedFig_EUD(loga,4,'r*',1);
        xlim([0 60]);
           
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_resp_m038'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_resp_m038.png...']);
        end
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        loga=CGnki.fLymanResponseFig_m_EUD(0.40);
        loga=CGnki.fComplicationObservedFig_EUD(loga,4,'r*',1);
        xlim([0 60]);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_resp_m040'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_resp_m040.png...']);
        end
        
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        loga=CGnki.fLymanResponseFig_m_EUD(0.42);
        CGnki.fComplicationObservedFig_EUD(loga,4,'r*',1);
        xlim([0 60]);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_resp_m042'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_resp_m042.png...']);
        end
        
        cur_fig=figure(6); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        loga=CGnki.fLymanResponseFig_m_EUD(0.44);
        CGnki.fComplicationObservedFig_EUD(loga,4,'r*',1);
        xlim([0 60]);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_resp_m044'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_resp_m044.png...']);
        end

        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanResponseFig_m_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_resp_m'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_resp_m.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanResponseFig_m_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_resp_m'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_resp_m.png...']);
        end
        
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanResponseFig_m_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_resp_m'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_resp_m.png...']);
        end
    
end

%maps - TD50 & a %toprint
if 0
          
    %disp(' '); disp(Regions{k});
          if ~is_comb_only, 
              
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGmsk.fLymanGridExactFig_TD50_a_EUD();
                
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_cntr_ll_td50_a'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_cntr_ll_td50_a.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGnki.fLymanGridExactFig_TD50_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_cntr_ll_td50_a'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_cntr_ll_td50_a.png...']);
        end
        
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanGridExactFig_TD50_a_EUD();

        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_cntr_ll_td50_a'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_cntr_ll_td50_a.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanGridExactFig_TD50_a_EUD();
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_cntr_ll_td50_a'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_cntr_ll_td50_a.png...']);
        end

          end
          
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_TD50_a_EUD();
      
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_cntr_ll_td50_a'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_cntr_ll_td50_a.png...']);
        end
        
        
        cur_fig=figure(6); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_TD50_a_EUD();hold on;
        xlim([-0.8, 0.2]);
        ylim([12,32]);
        plot(re_log10a(1),re_td50(1),'mx','MarkerSize',15,'LineWidth',2);
        h_rem=errorbar(re_log10a(1),re_td50(1),re_td50(1)-re_td50(2),re_td50(3)-re_td50(1),'m.-','LineWidth',2);
        h_rem_ex=errorbar_x(re_log10a(1),re_td50(1),re_log10a(1)-re_log10a(2),re_log10a(3)-re_log10a(1),'m.-');
        set(h_rem_ex,'LineWidth',2);
       
        plot(fe_log10a(1),fe_td50(1),'gx','MarkerSize',15,'LineWidth',2);
        h_fem=errorbar(fe_log10a(1),fe_td50(1),fe_td50(1)-fe_td50(2),fe_td50(3)-fe_td50(1),'g.-','LineWidth',2);
        h_fem_ex=errorbar_x(fe_log10a(1),fe_td50(1),fe_log10a(1)-fe_log10a(2),fe_log10a(3)-fe_log10a(1),'g.-');
        set(h_fem_ex,'LineWidth',2);
       
        set(gca,'box','on');
        h_lgnd=legend([h_fem h_rem],'Fixed Effects (95% CI)','Random Effects (95% CI)','Location','Best');
        set(h_lgnd,'FontSize',18);
      
      
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_td50_a_with_meta'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_td50_a_with_meta.png...']);
        end

end


%maps - m & a toprint
if 0
          
    %disp(' '); disp(Regions{k});
      if ~is_comb_only,        
        cur_fig=figure(1); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGmsk.fLymanGridExactFig_m_a_EUD();
        set(gca,'YLim',[0,1.2]);                
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'msk_lkb_cntr_ll_m_a'],'-png');
            disp(['Saving ',fig_loc,'msk_lkb_cntr_ll_m_a.png...']);
        end
        
          
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGnki.fLymanGridExactFig_m_a_EUD();
        set(gca,'YLim',[0,1.2]);                
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'nki_lkb_cntr_ll_m_a'],'-png');
            disp(['Saving ',fig_loc,'nki_lkb_cntr_ll_m_a.png...']);
        end
        
        
        cur_fig=figure(3); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGrtog.fLymanGridExactFig_m_a_EUD();
        set(gca,'YLim',[0,1.2]);                
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'rtog_lkb_cntr_ll_m_a'],'-png');
            disp(['Saving ',fig_loc,'rtog_lkb_cntr_ll_m_a.png...']);
        end
        
        
        cur_fig=figure(4); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGum.fLymanGridExactFig_m_a_EUD();
        set(gca,'YLim',[0,1.2]);                
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'umich_lkb_cntr_ll_m_a'],'-png');
            disp(['Saving ',fig_loc,'umich_lkb_cntr_ll_m_a.png...']);
        end
    
      end
        
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_m_a_EUD();
        set(gca,'YLim',[0,1.2]);                
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_lkb_cntr_ll_m_a'],'-png');
            disp(['Saving ',fig_loc,'comb_lkb_cntr_ll_m_a.png...']);
        end
    
        cur_fig=figure(6); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        CGcomb.fLymanGridExactFig_m_a_EUD();hold on;
        
         xlim([-0.8, 0.2]);
        ylim([0.1,0.6]);
    
        plot(re_log10a(1),re_m(1),'mx','MarkerSize',15,'LineWidth',2);
        h_rem=errorbar(re_log10a(1),re_m(1),re_m(1)-re_m(2),re_m(3)-re_m(1),'m.-','LineWidth',2);
        h_rem_ex=errorbar_x(re_log10a(1),re_m(1),re_log10a(1)-re_log10a(2),re_log10a(3)-re_log10a(1),'m.-');
        set(h_rem_ex,'LineWidth',2);
       
        plot(fe_log10a(1),fe_m(1),'gx','MarkerSize',15,'LineWidth',2);
        h_fem=errorbar(fe_log10a(1),fe_m(1),fe_m(1)-fe_m(2),fe_m(3)-fe_m(1),'g.-','LineWidth',2);
        h_fem_ex=errorbar_x(fe_log10a(1),fe_m(1),fe_log10a(1)-fe_log10a(2),fe_log10a(3)-fe_log10a(1),'g.-');
        set(h_fem_ex,'LineWidth',2);
       
        set(gca,'box','on');
        h_lgnd=legend([h_fem h_rem],'Fixed Effects (95% CI)','Random Effects (95% CI)','Location','Best');
        set(h_lgnd,'FontSize',18);

         if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_ll_m_a_with_meta'],'-png');
            disp(['Saving ',fig_loc,'comb_ll_m_a_with_meta.png...']);
        end
    
end


% maps - TD50 & a
% if 0
%     
% %     figure(1); clf reset; colormap(cm2);
% %     CGtmp.fLymanGridExactFig_TD50_a_EUD();
% %     title('TMP','FontSize',14);
% %     %set(gca,'YLim',[0.1,2.2]);
% %     xlabel(''); ylabel('');
% %     set(gca,'box','on');
%     
%     figure(1); clf reset; colormap(cm2);
%     CGmsk.fLymanGridExactFig_TD50_a_EUD();
%     title('MSK','FontSize',14);
%     %xlabel(''); ylabel('');
%     %set(gca,'box','on');
%     
%     figure(2); clf reset; colormap(cm2);
%     CGnki.fLymanGridExactFig_TD50_a_EUD();
%     title('NKI','FontSize',14);
%     %xlabel(''); ylabel('');
%     %set(gca,'box','on');
%     
%     figure(3); clf reset; colormap(cm2);
%     CGrtog.fLymanGridExactFig_TD50_a_EUD();
%     title('RTOG','FontSize',14);
%     %xlabel(''); ylabel('');
%     %set(gca,'box','on');
%     
%     figure(4); clf reset; colormap(cm2);
%     CGum.fLymanGridExactFig_TD50_a_EUD();
%     title('UMich','FontSize',14);
%     %xlabel(''); ylabel('');
%     %set(gca,'box','on');
%     
%     figure(5); clf reset; colormap(cm2);
%     CGcomb.fLymanGridExactFig_TD50_a_EUD();
%     title('COMB','FontSize',14);
%     %xlabel(''); ylabel('');
%     %set(gca,'box','on');
%     
% end

% maps - m & a
% if 0
%     
%     figure(1); clf reset; colormap(cm2);
%     CGmsk.fLymanGridExactFig_m_a_EUD();
%     title('MSK','FontSize',14);
%     set(gca,'YLim',[0.1,2.2]);
%     xlabel(''); ylabel('');
%     set(gca,'box','on');
%     
%     figure(2); clf reset; colormap(cm2);
%     CGnki.fLymanGridExactFig_m_a_EUD();
%     title('NKI','FontSize',14);
%     set(gca,'YLim',[0.1,2.2]);
%     xlabel(''); ylabel('');
%     set(gca,'box','on');
%     
%     figure(3); clf reset; colormap(cm2);
%     CGrtog.fLymanGridExactFig_m_a_EUD();
%     title('RTOG','FontSize',14);
%     set(gca,'YLim',[0.1,2.2]);
%     xlabel(''); ylabel('');
%     set(gca,'box','on');
%     
%     figure(4); clf reset; colormap(cm2);
%     CGum.fLymanGridExactFig_m_a_EUD();
%     title('UMich','FontSize',14);
%     set(gca,'YLim',[0.1,2.2]);
%     xlabel(''); ylabel('');
%     set(gca,'box','on');
%     
%     figure(5); clf reset; colormap(cm2);
%     CGcomb.fLymanGridExactFig_m_a_EUD();
%     title('COMB','FontSize',14);
%     set(gca,'YLim',[0.1,2.2]);
%     xlabel(''); ylabel('');
%     set(gca,'box','on');
%     
% end

% lyman curves - a
if 0

%         %%MSK
%         cur_fig=figure(1); clf reset; colormap(cm2);
%         set(gcf,'Position',ss_four2three);
%         h_msk=CGmsk.fLymanGridExactFig_a_loglikelihood('rs--',2);
%         h_lgnd=legend(h_msk,'MSK','Location','SouthWest');
%         ylim([-0.55 -0.2])
%         set(h_lgnd,'FontSize',18);
%         if do_print,
%             set(cur_fig,'Color','w');
%             export_fig(cur_fig,[fig_loc,'msk_lkb_llhds'],'-pdf');
%             disp(['Saving ',fig_loc,'msk_lkb_llhds.pdf...']);
%         end
%         
%         %%MSK + NKI
%         cur_fig=figure(2); clf reset; colormap(cm2);
%         set(gcf,'Position',ss_four2three);
%         h_msk=CGmsk.fLymanGridExactFig_a_loglikelihood('rs--',2);
%         h_nki=CGnki.fLymanGridExactFig_a_loglikelihood('bs--',2);
%         h_lgnd=legend([h_msk h_nki],'MSK','NKI','Location','SouthWest');
%         ylim([-0.55 -0.2])
%         set(h_lgnd,'FontSize',18);
%         if do_print,
%             set(cur_fig,'Color','w');
%             export_fig(cur_fig,[fig_loc,'msk_nki_lkb_llhds'],'-pdf');
%             disp(['Saving ',fig_loc,'msk_nki_lkb_llhds.pdf...']);
%         end
%         
%         %%MSK + NKI + UMich
%         cur_fig=figure(3); clf reset; colormap(cm2);
%         set(gcf,'Position',ss_four2three);
%         h_msk=CGmsk.fLymanGridExactFig_a_loglikelihood('rs--',2);
%         h_nki=CGnki.fLymanGridExactFig_a_loglikelihood('bs--',2);
%         h_um=CGum.fLymanGridExactFig_a_loglikelihood('ms--',2);
%         
%         
%         h_lgnd=legend([h_msk h_nki h_um],'MSK','NKI','UMich','Location','SouthWest');
%         ylim([-0.55 -0.2])
%         set(h_lgnd,'FontSize',18);
%         if do_print,
%             set(cur_fig,'Color','w');
%             export_fig(cur_fig,[fig_loc,'msk_nki_um_lkb_llhds'],'-pdf');
%             disp(['Saving ',fig_loc,'msk_nki_um_lkb_llhds.pdf...']);
%         end
%         
%         
%            %%MSK + NKI + UMich + RTOG
%         cur_fig=figure(4); clf reset; colormap(cm2);
%         set(gcf,'Position',ss_four2three);
%         h_msk=CGmsk.fLymanGridExactFig_a_loglikelihood('rs--',2);
%         h_nki=CGnki.fLymanGridExactFig_a_loglikelihood('bs--',2);
%         h_um=CGum.fLymanGridExactFig_a_loglikelihood('ms--',2);
%         h_rtog=CGrtog.fLymanGridExactFig_a_loglikelihood('gs--',2);
%         h_lgnd=legend([h_msk h_nki h_um h_rtog],'MSK','NKI','UMich','RTOG','Location','SouthWest');
%         ylim([-0.55 -0.2])
%         set(h_lgnd,'FontSize',18);
%         if do_print,
%             set(cur_fig,'Color','w');
%             export_fig(cur_fig,[fig_loc,'msk_nki_um_rtog_lkb_llhds'],'-pdf');
%             disp(['Saving ',fig_loc,'msk_nki_um_rtog_lkb_llhds.pdf...']);
%         end
        
        %%MSK + NKI + UMich + RTOG + COMB
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        disp('******');
        disp([]);
        disp('MSK');
        h_msk=CGmsk.fLymanGridExactFig_a_loglikelihood('r-',2);
        disp('NKI');
        h_nki=CGnki.fLymanGridExactFig_a_loglikelihood('b-',2);
        disp('UM');
        h_um=CGum.fLymanGridExactFig_a_loglikelihood('m-',2);
        disp('RTOG');
        h_rtog=CGrtog.fLymanGridExactFig_a_loglikelihood('g-',2);
        disp('COMB');
        h_comb=CGcomb.fLymanGridExactFig_a_loglikelihood('k-',2);
        
      

    
        
%         xlabel(''); ylabel('');
%          set(gca,'box','on');
% %         % set log10(n) at the top
%           a1 = gca;
%           a2 = copyobj(a1,gcf);
%           set(a2,'Color','none');
%           set(a2,'XAxisLocation','top');
%           set(a2,'XTickLabel',num2str(CGcomb.mLymanN(end:-2:1)));
%          
%           set(get(a1,'XLabel'),'String','log_1_0(a)');
%           set(get(a1,'XLabel'),'FontSize',20);
%           
        set(gca,'FontSize',18);
        ylabel('Log likelihood per degree of freedom','FontSize',20);
        xlabel('log_{10}(a)','FontSize',20);
              
        ylim([-0.55 -0.2]);
              
        h_lgnd=legend([h_msk h_nki h_um h_rtog h_comb],'MSK','NKI','UMich','RTOG','COMB','Location','SouthWest');
        set(h_lgnd,'FontSize',18);
        %
        %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'lkb_ll_log10a'],'-png');
            disp(['Saving ',fig_loc,'lkb_ll_log10a.png...']);
        end


end
% lyman curves - n
if 0


        %%MSK + NKI + UMich + RTOG + COMB
        cur_fig=figure(5); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        disp('******');
        disp([]);
        disp('MSK');
        h_msk=CGmsk.fLymanGridExactFig_n_loglikelihood('r-',2);
        disp('NKI');
        h_nki=CGnki.fLymanGridExactFig_n_loglikelihood('b-',2);
        disp('UM');
        h_um=CGum.fLymanGridExactFig_n_loglikelihood('m-',2);
        disp('RTOG');
        h_rtog=CGrtog.fLymanGridExactFig_n_loglikelihood('g-',2);
        disp('COMB');
        h_comb=CGcomb.fLymanGridExactFig_n_loglikelihood('k-',2);
  
        set(gca,'FontSize',20);
        ylabel('Log likelihood per degree of freedom','FontSize',24);
        xlabel('log_{10}(n)','FontSize',24);
              
        ylim([-0.55 -0.2]);
              
        h_lgnd=legend([h_msk h_nki h_um h_rtog h_comb],'MSK','NKI','UMich','RTOG','COMB','Location','SouthWest');
        set(h_lgnd,'FontSize',18);
        set(h_lgnd,'Location','SouthEast')
        %
        %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'lkb_ll_log10n'],'-png');
            disp(['Saving ',fig_loc,'lkb_ll_log10n.png...']);
        end


end
% m
if 0
%        disp(' '); disp(Regions{k});
        cur_fig=figure(6); clf reset; colormap(cm2);
        set(cur_fig,'Position',ss_four2three);
        disp('********');
        disp('MSK');
        h_msk=CGmsk.fLymanGridExactFig_m_loglikelihood('r',2);
        disp([]);
        disp('NKI');
        h_nki=CGnki.fLymanGridExactFig_m_loglikelihood('b',2);
        
        disp([]);
        disp('RTOG');
        h_rtog=CGrtog.fLymanGridExactFig_m_loglikelihood('g',2);
        
        disp([]);
        disp('UMich');
        h_um=CGum.fLymanGridExactFig_m_loglikelihood('m',2);
        
        disp([]);
        disp('Comb');
        CGcomb.fLymanGridExactFig_m_loglikelihood('k--',2);
        %h_comb=CGcomb.fLymanGridExactFig_m_loglikelihoodAtLoga_EUD(-0.2,'ks--',2);
        %        CGcomb.fLymanGridExactFig_m_loglikelihoodAtLoga_EUD(0,'ks--',2);
        %CGcomb.fLymanGridExactFig_m_loglikelihood(0,'ks--',2);
        
        %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
        
        set(gca,'FontSize',20);
        ylabel('Log likelihood per degree of freedom','FontSize',24);
        xlabel('m','FontSize',24);
        ylim([-0.5 -0.2]);
        xlim([0.1 1.2]);
        %set(h_lgnd,'FontSize',15);
        %if do_print, print_fig(gcf,fig_loc,'lkb_ll_m','png');end;
        
        if do_print
             set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'lkb_ll_m'],'-png');
            disp(['Saving ',fig_loc,'lkb_ll_m.png...']);
        end
        
%         xlabel(''); ylabel('');
%         set(gca,'box','on');
%         % set log10(n) at the top
%         a1 = gca;
%         a2 = copyobj(a1,gcf);
%         set(a2,'Color','none');
%         set(a2,'XAxisLocation','top');
%         set(a2,'XTickLabel',num2str(CGcomb.mLymanN(end:-2:1)));

end

% TD50
if 0%        disp(' '); disp(Regions{k});
        cur_fig=figure(2); clf reset; colormap(cm2);
        set(cur_fig,'Position',ss_four2three);
        disp('*******');
        disp('MSK');
         h_msk=CGmsk.fLymanGridExactFig_TD50_loglikelihood('r',2);
         disp('NKI');
         h_nki=CGnki.fLymanGridExactFig_TD50_loglikelihood('b',2);
         disp('RTOG');
         h_rtog=CGrtog.fLymanGridExactFig_TD50_loglikelihood('g',2);
         disp('UMich');
         h_um=CGum.fLymanGridExactFig_TD50_loglikelihood('m',2);
         
%        CGtmp.fLymanGridExactFig_TD50_loglikelihood('ms--',2);
        disp('Comb');
        CGcomb.fLymanGridExactFig_TD50_loglikelihood('k',2);
        
        
        %h_comb=CGcomb.fLymanGridExactFig_TD50_loglikelihoodAtLoga_EUD(-0.2,'ks--',2);
        
       
           set(gca,'FontSize',20);
        ylabel('Log likelihood per degree of freedom','FontSize',24);
        xlabel('TD_{50}','FontSize',24);
        
        ylim([-0.5 -0.2]);
        xlim([0 72]);
        %set(h_lgnd,'FontSize',15);
        %if do_print, print_fig(gcf,fig_loc,'lkb_ll_m','png');end;
        
        if do_print
             set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'lkb_ll_td50'],'-png');
            disp(['Saving ',fig_loc,'lkb_ll_td50.png...']);
        end

%        disp([]);
%        disp('******');
%        disp('Comb MLD');
%        CGcomb.fLymanGridExactFig_TD50_loglikelihoodAtLoga_EUD(0,'ks--',2);



end

% response curve
if 0
    disp(' '); disp(['Whole']);
    figure(10); clf reset; colormap(cm2);
    loga = CGtmp.fLymanGridResponseExactFig_a_EUD('loga','r',1);
    CGtmp.fComplicationObservedFig_EUD(loga,4,'r*',1);
    if do_print, print_fig(gcf,fig_loc,'tmp_lkb_reponse','png');end;
    xlabel(''); ylabel('');
    set(gca,'box','on');
   
    figure(1); clf reset; colormap(cm2);
    loga = CGmsk.fLymanGridResponseExactFig_a_EUD('loga','r',1);
    CGmsk.fComplicationObservedFig_EUD(loga,4,'r*',1);
    if do_print, print_fig(gcf,fig_loc,'msk_lkb_reponse','png');end;
    xlabel(''); ylabel('');
    set(gca,'box','on');
    
    figure(2); clf reset; colormap(cm2);
    loga = CGnki.fLymanGridResponseExactFig_a_EUD('loga','b',1);
    CGnki.fComplicationObservedFig_EUD(loga,4,'b*',1);
    if do_print, print_fig(gcf,fig_loc,'nki_lkb_reponse','png');end;
    xlabel(''); ylabel('');
    set(gca,'box','on');
      
    figure(3); clf reset; colormap(cm2);
    loga = CGrtog.fLymanGridResponseExactFig_a_EUD('loga','g',1);
    CGrtog.fComplicationObservedFig_EUD(loga,4,'g*',1);
    if do_print, print_fig(gcf,fig_loc,'rtog_lkb_reponse','png');end;
    xlabel(''); ylabel('');
    set(gca,'box','on');
    
    cur_fig=figure(4); clf reset; colormap(cm2);
    loga = CGcomb.fLymanGridResponseExactFig_a_EUD('loga','k',1);
    CGcomb.fComplicationObservedFig_EUD(loga,4,'k*',1);
    %if do_print, print_fig(gcf,fig_loc,'comb_lkb_reponse','png');end;
    set(gca,'YLim',[0,0.8])
    set(gca,'box','on');
    if do_print
         set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'comb_lkb_response'],'-png');
        disp(['Saving ',fig_loc,'comb_lkb_response.png...']);
    end

    
end

% data sets consistency, heterogeniety and inconsistency
if 0
    %% MSK + NKI
    disp(' '); 
%     disp('MSK + NKI: Consistency LKB');
%     [Q,Qp, I2,I2up,I2down] = classOutcomeAnalysis.fLKBDataConsistency_EUD([CGmsk; CGnki]');
%     disp(['Var','   Q  ',' Qp ',' I2 ',' I2up ',' I2down ']);
%     %disp([Q,Qp,I2,I2up,I2down]);% Q (TD50, m, a)
%     
%      disp(['TD50 & ',...
%         num2str(Q(1),4),' & ',...
%         num2str(Qp(1),3),' & ',...
%         num2str(I2(1),3),' & ',...
%         num2str(I2up(1),3),' & ',...
%         num2str(I2down(1),3)]);
%     disp(['m & ',...
%         num2str(Q(2),4),' & ',...
%         num2str(Qp(2),3),' & ',...
%         num2str(I2(2),3),' & ',...
%         num2str(I2up(2),3),' & ',...
%         num2str(I2down(2),3)]);
%     disp(['a & ',...
%         num2str(Q(3),4),' & ',...
%         num2str(Qp(3),3),' & ',...
%         num2str(I2(3),3),' & ',...
%         num2str(I2up(3),3),' & ',...
%         num2str(I2down(3),3)]);
%    
%     %% MSK + NKI + RTOG
%     disp(' '); 
%     disp('MSK + NKI + RTOG: Consistency LKB');
%     [Q,Qp, I2,I2up,I2down] = classOutcomeAnalysis.fLKBDataConsistency_EUD([CGmsk; CGnki;CGrtog]');
%     disp(['Var','   Q  ',' Qp ',' I2 ',' I2up ',' I2down ']);
%     %disp([Q,Qp,I2,I2up,I2down]);% Q (TD50, m, a)
%     
%      disp(['TD50 & ',...
%         num2str(Q(1),4),' & ',...
%         num2str(Qp(1),3),' & ',...
%         num2str(I2(1),3),' & ',...
%         num2str(I2up(1),3),' & ',...
%         num2str(I2down(1),3)]);
%     disp(['m & ',...
%         num2str(Q(2),4),' & ',...
%         num2str(Qp(2),3),' & ',...
%         num2str(I2(2),3),' & ',...
%         num2str(I2up(2),3),' & ',...
%         num2str(I2down(2),3)]);
%     disp(['a & ',...
%         num2str(Q(3),4),' & ',...
%         num2str(Qp(3),3),' & ',...
%         num2str(I2(3),3),' & ',...
%         num2str(I2up(3),3),' & ',...
%         num2str(I2down(3),3)]);
% 
%      %% MSK + NKI + UMich
%     disp(' '); 
%     disp('MSK + NKI + UMich: Consistency LKB');
%     [Q,Qp, I2,I2up,I2down] = classOutcomeAnalysis.fLKBDataConsistency_EUD([CGmsk; CGnki;CGum]');
%     disp(['Var','   Q  ',' Qp ',' I2 ',' I2up ',' I2down ']);
%     %disp([Q,Qp,I2,I2up,I2down]);% Q (TD50, m, a)
%     
%      disp(['TD50 & ',...
%         num2str(Q(1),4),' & ',...
%         num2str(Qp(1),3),' & ',...
%         num2str(I2(1),3),' & ',...
%         num2str(I2up(1),3),' & ',...
%         num2str(I2down(1),3)]);
%     disp(['m & ',...
%         num2str(Q(2),4),' & ',...
%         num2str(Qp(2),3),' & ',...
%         num2str(I2(2),3),' & ',...
%         num2str(I2up(2),3),' & ',...
%         num2str(I2down(2),3)]);
%     disp(['a & ',...
%         num2str(Q(3),4),' & ',...
%         num2str(Qp(3),3),' & ',...
%         num2str(I2(3),3),' & ',...
%         num2str(I2up(3),3),' & ',...
%         num2str(I2down(3),3)]);
     %% MSK + NKI + UMich + RTOG
    disp(' '); 
    disp('MSK + NKI + RTOG + UMich: Consistency LKB');
    [Q,Qp, I2,I2up,I2down] = classOutcomeAnalysis.fLKBDataConsistency_EUD([CGmsk; CGnki;CGum;CGrtog]');
    disp(['Var','   Q  ',' Qp ',' I2 ',' I2up ',' I2down ']);
    %disp([Q,Qp,I2,I2up,I2down]);% Q (TD50, m, a)
    
     disp(['TD50 & ',...
        num2str(Q(1),4),' & ',...
        num2str(Qp(1),3),' & ',...
        num2str(I2(1),3),' & ',...
        num2str(I2up(1),3),' & ',...
        num2str(I2down(1),3)]);
    disp(['m & ',...
        num2str(Q(2),4),' & ',...
        num2str(Qp(2),3),' & ',...
        num2str(I2(2),3),' & ',...
        num2str(I2up(2),3),' & ',...
        num2str(I2down(2),3)]);
    disp(['a & ',...
        num2str(Q(3),4),' & ',...
        num2str(Qp(3),3),' & ',...
        num2str(I2(3),3),' & ',...
        num2str(I2up(3),3),' & ',...
        num2str(I2down(3),3)]);
end

% Logistic Regression
% data sets consistency
if 0
    %for k = 1:numreg
    disp(' '); 
    disp('Consistency Logistic Regression');
    [Q,Qp, I2] = classOutcomeAnalysis.fLogisticRegressionDataConsistency_EUD([CGmsk; CGnki;CGum;CGrtog]');
    disp(['  Q ','  Qp ','  I2  ']);
    disp(['TD50: ',...
        num2str(Q(1),3),' ',...
        num2str(Qp(1),3),' ',...
        num2str(I2(1),3)]);
    disp(['m: ',...
        num2str(Q(2),3),' ',...
        num2str(Qp(2),3),' ',...
        num2str(I2(2),3)]);
    disp(['a: ',...
        num2str(Q(3),3),' ',...
        num2str(Qp(3),3),' ',...
        num2str(I2(3),3)]);
        
end
% p-values from matlab function
if 0
    
    figure(60); clf reset; colormap(cm2);
    set(gcf,'Position',ss_four2three);
    h_msk=CGmsk.fLogisticRegressionPvalueExactFig_a_EUD('rs--',2);
    h_nki=CGnki.fLogisticRegressionPvalueExactFig_a_EUD('bs--',2);
    h_rtog=CGrtog.fLogisticRegressionPvalueExactFig_a_EUD('gs--',2);
    h_um=CGum.fLogisticRegressionPvalueExactFig_a_EUD('ms--',2);    
    h_comb=CGcomb.fLogisticRegressionPvalueExactFig_a_EUD('ks--',2);
    h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
    set(h_lgnd,'FontSize',15);
    if do_print, print_fig(gcf,fig_loc,'logreg_pvals','png');end;
    xlabel(''); ylabel('');
    set(gca,'box','on');
    
end

% lyman analysis p-values
if 0

    cur_fig=figure(71); clf reset; colormap(cm2);
    set(gcf,'Position',ss_four2three);
    h_msk=CGmsk.fLKBPvalueFig_a_EUD('rs--',2);
    h_lgnd=legend(h_msk,'MSK','Location','SouthEast');
    ylim([10^(-10) 1]);
    %h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    %h_um=CGum.fLKBPvalueFig_a_EUD('ms--',2);    
    %h_rtog=CGrtog.fLKBPvalueFig_a_EUD('gs--',2);
    %h_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
    %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
    set(h_lgnd,'FontSize',18);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'msk_lkb_pvals'],'-pdf');
        disp(['Saving ',fig_loc,'msk_lkb_pvals.pdf...']);
    end
    
    cur_fig=figure(72); clf reset; colormap(cm2);
    set(gcf,'Position',ss_four2three);
    h_msk=CGmsk.fLKBPvalueFig_a_EUD('rs--',2);
    h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    h_lgnd=legend([h_msk h_nki],'MSK','NKI','Location','SouthEast');
    ylim([10^(-10) 1]);
    %h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    %h_um=CGum.fLKBPvalueFig_a_EUD('ms--',2);    
    %h_rtog=CGrtog.fLKBPvalueFig_a_EUD('gs--',2);
    %h_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
    %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
    set(h_lgnd,'FontSize',18);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'msk_nki_lkb_pvals'],'-pdf');
        disp(['Saving ',fig_loc,'msk_nki_lkb_pvals.pdf...']);
    end
    
    %% MSK NKI UMich
    cur_fig=figure(73); clf reset; colormap(cm2);
    set(gcf,'Position',ss_four2three);
    h_msk=CGmsk.fLKBPvalueFig_a_EUD('rs--',2);
    h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    h_um=CGum.fLKBPvalueFig_a_EUD('ms--',2);    
    h_lgnd=legend([h_msk h_nki h_um],'MSK','NKI','UMich','Location','SouthEast');
    ylim([10^(-10) 1]);
    
    %h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    
    %h_rtog=CGrtog.fLKBPvalueFig_a_EUD('gs--',2);
    %h_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
    %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
    set(h_lgnd,'FontSize',18);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'msk_nki_um_lkb_pvals'],'-pdf');
        disp(['Saving ',fig_loc,'msk_nki_um_lkb_pvals.pdf...']);
    end
    
    %% MSK NKI UMich RTOG
    cur_fig=figure(74); clf reset; colormap(cm2);
    set(gcf,'Position',ss_four2three);
    h_msk=CGmsk.fLKBPvalueFig_a_EUD('rs--',2);
    h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    h_um=CGum.fLKBPvalueFig_a_EUD('ms--',2);    
    h_rtog=CGrtog.fLKBPvalueFig_a_EUD('gs--',2);
    h_lgnd=legend([h_msk h_nki h_um h_rtog],'MSK','NKI','UMich','RTOG','Location','SouthEast');
    ylim([10^(-10) 1]);

    %h_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
    %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
    set(h_lgnd,'FontSize',18);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'msk_nki_um_rtog_lkb_pvals'],'-pdf');
        disp(['Saving ',fig_loc,'msk_nki_um_rtog_lkb_pvals.pdf...']);
    end
%     
  %% MSK NKI UMich RTOG
    cur_fig=figure(75); clf reset; colormap(cm2);
    set(gcf,'Position',ss_four2three);
    h_msk=CGmsk.fLKBPvalueFig_a_EUD('rs--',2);
    h_nki=CGnki.fLKBPvalueFig_a_EUD('bs--',2);
    h_um=CGum.fLKBPvalueFig_a_EUD('ms--',2);    
    h_rtog=CGrtog.fLKBPvalueFig_a_EUD('gs--',2);
    h_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
    h_lgnd=legend([h_msk h_nki h_um h_rtog h_comb],'MSK','NKI','UMich','RTOG','COMB','Location','SouthEast');
    ylim([10^(-10) 1]);

    %
    %h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
    set(h_lgnd,'FontSize',18);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'comb_lkb_pvals'],'-pdf');
        disp(['Saving ',fig_loc,'comb_lkb_pvals.pdf...']);
    end
end


% lyman analysis p-values
if 0
    ctr=1;
    for i=length(LymanN):-1:1,
        cur_loga=-LymanN(i);
        cur_fig=figure(i); clf reset; colormap(cm2);
        set(gcf,'Position',2.*ss_four2three);

        subplot(2,4,[1:2]);
        h_llhd_comb=CGcomb.fLymanGridExactFig_a_loglikelihood('ks--',2);
        line([cur_loga cur_loga],ylim,'LineWidth',3,'Color','r');
        
        subplot(2,4,[3:4]);
        h_pval_comb=CGcomb.fLKBPvalueFig_a_EUD('ks--',2);
        line([cur_loga cur_loga],ylim,'LineWidth',3,'Color','r');
    
  
        subplot(2,4,[6:7])
        [loga,pval] = CGcomb.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'r',1);
        CGcomb.fComplicationObservedFig_EUD(loga,4,'r*',2);
        loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
        str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
        text(0.65,0.85,str,'FontSize',24,'Units','normalized');
        %text(2,0.6,str,'FontSize',18);
        title('MSK + NKI + RTOG + UMich','FontSize',24);    
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'lkb_comb2-',num2str(ctr)],'-pdf');
        disp(['Saving ',fig_loc,'lkb_comb2-',num2str(ctr),'.pdf...']);
    end
    ctr=ctr+1;
    end
end    



% maps - b0 & b1
if 0
    for k = 1:numreg
        disp(' '); disp(Regions{k});
        
        figure(1); clf reset; colormap(cm2);
        CGmsk(k).fLogisticRegressionGridExactFig_b0_b1_EUD();
        set(gca,'YLim',[-0.6,0.9]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
        figure(2); clf reset; colormap(cm2);
        CGnki(k).fLogisticRegressionGridExactFig_b0_b1_EUD();
        set(gca,'YLim',[-0.6,0.9]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
        figure(3); clf reset; colormap(cm2);
        CGcomb(k).fLogisticRegressionGridExactFig_b0_b1_EUD();
        set(gca,'YLim',[-0.6,0.9]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
    end
end
    
% maps - b0 & a
if 0
    for k = 1:numreg
        disp(' '); disp(Regions{k});
        figure(1); clf reset; colormap(cm2);
        CGmsk(k).fLogisticRegressionGridExactFig_b0_a_EUD();
        set(gca,'YLim',[-10,4]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        figure(2); clf reset; colormap(cm2);
        CGnki(k).fLogisticRegressionGridExactFig_b0_a_EUD();
        set(gca,'YLim',[-10,4]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        figure(3); clf reset; colormap(cm2);
        CGcomb(k).fLogisticRegressionGridExactFig_b0_a_EUD();
        set(gca,'YLim',[-10,4]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
    end
end

% maps - b1 & a
if 0
    for k = 1:numreg
        disp(' '); disp(Regions{k});
        figure(1); clf reset; colormap(cm2);
        CGmsk(k).fLogisticRegressionGridExactFig_b1_a_EUD();
        set(gca,'YLim',[-0.4,0.7]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        figure(2); clf reset; colormap(cm2);
        CGnki(k).fLogisticRegressionGridExactFig_b1_a_EUD();
        set(gca,'YLim',[-0.4,0.7]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
        figure(3); clf reset; colormap(cm2);
        CGcomb(k).fLogisticRegressionGridExactFig_b1_a_EUD();
        set(gca,'YLim',[-0.4,0.7]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
    end
end

% curves - a
if 0
        %disp(' '); disp(Regions{k});
        figure(7); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        h_msk=CGmsk.fLogisticRegressionGridExactFig_a_loglikelhood_EUD('rs--',2);
        h_nki=CGnki.fLogisticRegressionGridExactFig_a_loglikelhood_EUD('bs--',2);
        h_um=CGum.fLogisticRegressionGridExactFig_a_loglikelhood_EUD('ms--',2);
        h_rtog=CGrtog.fLogisticRegressionGridExactFig_a_loglikelhood_EUD('gs--',2);
        h_comb=CGcomb.fLogisticRegressionGridExactFig_a_loglikelhood_EUD('ks--',2);
        h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
        set(h_lgnd,'FontSize',15);
        if do_print, print_fig(gcf,fig_loc,'logreg_ll_loga','png');end;
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
        
end

% curves - b0
if 0
    for k = 1:numreg
        disp(' '); disp(Regions{k});
        figure(2); clf reset; colormap(cm2);
        CGmsk(k).fLogisticRegressionGridExactFig_b0_loglikelihood_EUD('r--',2);
        CGnki(k).fLogisticRegressionGridExactFig_b0_loglikelihood_EUD('b--',2);
        CGcomb(k).fLogisticRegressionGridExactFig_b0_loglikelihood_EUD('k--',2);
        set(gca,'XLim',[-10,10]); set(gca,'YLim',[-0.5,-0.3]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
    end
end

% curves - b1
if 0
    for k = 1:numreg
        disp(' '); disp(Regions{k});
        figure(3); clf reset; colormap(cm2);
        CGmsk(k).fLogisticRegressionGridExactFig_b1_loglikelihood_EUD('r--',2);
        CGnki(k).fLogisticRegressionGridExactFig_b1_loglikelihood_EUD('b--',2);
        CGcomb(k).fLogisticRegressionGridExactFig_b1_loglikelihood_EUD('k--',2);
        set(gca,'XLim',[-0.4,0.8]); set(gca,'YLim',[-0.5,-0.3]);
        xlabel(''); ylabel('');
        set(gca,'box','on');
    end
end

% response curve
if 0

        figure(1); clf reset; colormap(cm2);
        loga = CGmsk.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','r',1);
        CGmsk.fComplicationObservedFig_EUD(loga,4,'r*',1);
        if do_print, print_fig(gcf,fig_loc,'msk_logreg_reponse','png');end;
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
        figure(2); clf reset; colormap(cm2);
        loga = CGnki.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','b',1);
        CGnki.fComplicationObservedFig_EUD(loga,4,'b*',1);
        if do_print, print_fig(gcf,fig_loc,'nki_logreg_reponse','png');end;
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
        figure(3); clf reset; colormap(cm2);
        
        loga = CGrtog.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','g',1);
        CGrtog.fComplicationObservedFig_EUD(loga,4,'g*',1);
        if do_print, print_fig(gcf,fig_loc,'rtog_logreg_reponse','png');end;
        xlabel(''); ylabel('');
        set(gca,'box','on');
       
       
        
        figure(4); clf reset; colormap(cm2);
        loga = CGcomb.fLogisticRegressionRespondingCurveExactFig_a_EUD('loga','k',1);
        CGcomb.fComplicationObservedFig_EUD(loga,4,'k*',1);
        if do_print, print_fig(gcf,fig_loc,'comb_logreg_reponse','png');end;
        xlabel(''); ylabel('');
        set(gca,'box','on');
 
end
% 
%# Logistic Regression Response function for all loga
if 0
        for i=1:length(LymanN),
            cur_loga=-LymanN(i);
            figure(i);clf reset;colormap(cm2);
            set(gcf,'Position',ss_four2three);
            hold off;
            [loga,pval] = CGrtog.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'g',1);
            CGrtog.fComplicationObservedFig_EUD(loga,4,'g*',1);
            if loga~=cur_loga,
                disp(['loga ~= cur_loga!']);
            end
            loga_str = ['RTOG',10,'Log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,3));
            %text(2,0.6,str,'FontSize',16);
            text(0.8,0.85,str,'FontSize',16,'Units','normalized');
            
            fig_name = 'rtog_logreg_response_%s';
            fig_str = sprintf(fig_name,num2str(cur_loga));
            fig_str = strrep(fig_str,'.','_');
            fig_str = strrep(fig_str,'-','m');
            if do_print, print_fig(gcf,fig_loc,fig_str,'png');end;
            
            clf reset;colormap(cm2);
            [loga,pval] = CGmsk.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'r',1);
            CGmsk.fComplicationObservedFig_EUD(loga,4,'r*',1);
            loga_str = ['MSK',10,'Log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,3));
        %text(2,0.6,str,'FontSize',16);
            text(0.8,0.85,str,'FontSize',16,'Units','normalized');
            fig_name = 'msk_logreg_response_%s';
            fig_str = sprintf(fig_name,num2str(cur_loga));
            fig_str = strrep(fig_str,'.','_');
            fig_str = strrep(fig_str,'-','m');
            if do_print, print_fig(gcf,fig_loc,fig_str,'png');end;
            
            clf reset;colormap(cm2);
            [loga,pval] = CGnki.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'b',1);
            CGnki.fComplicationObservedFig_EUD(loga,4,'b*',1);
            loga_str = ['NKI',10,'Log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,3));            
            %text(2,0.6,str,'FontSize',16);
            text(0.8,0.85,str,'FontSize',16,'Units','normalized');
            fig_name = 'nki_logreg_response_%s';
            fig_str = sprintf(fig_name,num2str(cur_loga));
            fig_str = strrep(fig_str,'.','_');
            fig_str = strrep(fig_str,'-','m');
            if do_print, print_fig(gcf,fig_loc,fig_str,'png');end;
            
            clf reset;colormap(cm2);
            [loga,pval] = CGum.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'m',1);
            CGum.fComplicationObservedFig_EUD(loga,4,'m*',1);
            loga_str = ['UMich',10,'Log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,3));
            %text(2,0.6,str,'FontSize',16);
            text(0.8,0.85,str,'FontSize',16,'Units','normalized');
            fig_name = 'um_logreg_response_%s';
            fig_str = sprintf(fig_name,num2str(cur_loga));
            fig_str = strrep(fig_str,'.','_');
            fig_str = strrep(fig_str,'-','m');
            if do_print, print_fig(gcf,fig_loc,fig_str,'png');end;
    
            
            clf reset;colormap(cm2);
            [loga,pval] = CGcomb.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'k',1);
            CGcomb.fComplicationObservedFig_EUD(loga,4,'k*',1);
            loga_str = ['Comb',10,'Log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,3));
            %text(2,0.6,str,'FontSize',16);
            text(0.8,0.85,str,'FontSize',16,'Units','normalized');
            if do_print,
                %print_fig(gcf,fig_loc,fig_str,'png');end;
            set(gcf,'Color','w');
            export_fig(gcf,[fig_loc,'comb_logreg_response-',num2str(i)],'-pdf');
            end
        


        end
     
        
end


%# Logistic Regression Response function for all loga
% 4x4 with subplots
if 0
        for i=1:length(LymanN),
            cur_loga=-LymanN(i);
            cur_fig=figure(i);clf reset;colormap(cm2);
            set(cur_fig,'Position',[0 0 screen_size(3) screen_size(4)]);
            hold off;
            subplot(2,2,2)
            [loga,pval] = CGrtog.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'g',1);
            CGrtog.fComplicationObservedFig_EUD(loga,4,'g*',1);
            if loga~=cur_loga,
                disp(['loga ~= cur_loga!']);
            end
            loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
            text(0.75,0.85,str,'FontSize',16,'Units','normalized');
            %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
            %text(2,0.6,str,'FontSize',18);
            title('RTOG','FontSize',24);
            
            subplot(2,2,1)
            [loga,pval] = CGmsk.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'r',1);
            CGmsk.fComplicationObservedFig_EUD(loga,4,'r*',1);
            loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
            text(0.75,0.85,str,'FontSize',16,'Units','normalized');
            %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
            %text(2,0.6,str,'FontSize',18);
            title('MSK','FontSize',24);
            
            subplot(2,2,4)
            [loga,pval] = CGnki.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'b',1);
            CGnki.fComplicationObservedFig_EUD(loga,4,'b*',1);
            loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));            
            text(0.75,0.85,str,'FontSize',16,'Units','normalized');
            %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
            %text(2,0.6,str,'FontSize',18);
            title('NKI','FontSize',24);
            
            subplot(2,2,3)
            [loga,pval] = CGum.fLogisticRegressionRespondingCurveExactFig_a_EUD(cur_loga,'m',1);
            CGum.fComplicationObservedFig_EUD(loga,4,'m*',1);
            loga_str = ['log_1_0(a) = %s',10,'p-val: %s'];
            str = sprintf(loga_str,num2str(cur_loga),num2str(pval,2));
            text(0.75,0.85,str,'FontSize',16,'Units','normalized');
            %text(0.05,0.85,str,'FontSize',16,'Units','normalized');
            %text(2,0.6,str,'FontSize',18);
            title('UMich','FontSize',24);
            
            if do_print,
                set(cur_fig,'Color','w');
                export_fig(cur_fig,[fig_loc,'logreg_responses_frange',...
                    '-',num2str(i)],'-pdf');
                disp(['Saving ',fig_loc,'logreg_responses_frange',...
                    '-',num2str(i),'.pdf...']);
            end;
        end;
end;


%# LKB Response functions for all loga
if 0
    
    for i=1:length(LymanN),
        cur_loga=LymanN(i);
        figure(i); clf reset; colormap(cm2);
        loga = CGrtog.fLymanGridResponseExactFig_a_EUD(cur_loga,'r',1);
        CGrtog.fComplicationObservedFig_EUD(loga,4,'r*',1);
        
        loga_str = 'Log_1_0(a) = %s';
        if loga~=cur_loga,
            disp(['loga ~= cur_loga!']);
        end
        str = sprintf(loga_str,num2str(cur_loga));
        text(2,0.4,str,'FontSize',15);
            
        fig_name = 'rtog_lkb_response_%s';
        fig_str = sprintf(fig_name,num2str(cur_loga));
        fig_str = strrep(fig_str,'.','_');
        fig_str = strrep(fig_str,'-','m');
        if do_print, print_fig(gcf,fig_loc,fig_str,'png');end;
        
        
        
    end;
    
end

% log-likehood for mean doses

if 0
    cur_loga = 0;
    CGmsk.fLymanGridResponseExactFig_a_EUD(cur_loga,'r',1);
    CGnki.fLymanGridResponseExactFig_a_EUD(cur_loga,'r',1);
    CGrtog.fLymanGridResponseExactFig_a_EUD(cur_loga,'r',1);
    CGum.fLymanGridResponseExactFig_a_EUD(cur_loga,'r',1);    
    CGcomb.fLymanGridResponseExactFig_a_EUD(cur_loga,'r',1);    

end



if 0
        %disp(' '); disp(Regions{k});
        figure(80); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        h_msk=CGmsk.fSpearmanCorrelationCoefficient_a_EUD('rs--',2);
        h_nki=CGnki.fSpearmanCorrelationCoefficient_a_EUD('bs--',2);
        h_um=CGum.fSpearmanCorrelationCoefficient_a_EUD('ms--',2);
        h_rtog=CGrtog.fSpearmanCorrelationCoefficient_a_EUD('gs--',2);
        h_comb=CGcomb.fSpearmanCorrelationCoefficient_a_EUD('ks--',2);
        h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
        set(h_lgnd,'FontSize',15);
        if do_print, print_fig(gcf,fig_loc,'spearman_rho','png');end;
        xlabel(''); ylabel('');
        set(gca,'box','on');
        
            
    %    figure(2); clf reset; colormap(cm2);
    
    %    CGmsk.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga','rs--',2);
    %    CGnki.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga','bs--',2);
    %    CGrtog.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga','gs--',2);
    %    CGcomb.fLogisticRegressionLikelyhoodExactFig_a_EUD('loga','ks--',2);
    
    %    xlabel(''); ylabel('');
    %    set(gca,'box','on');
        
end


if 0
        figure(90); clf reset; colormap(cm2);
        set(gcf,'Position',ss_four2three);
        h_msk=CGmsk.fWilcoxonRankSum_a_EUD('rs--',2);
        h_nki=CGnki.fWilcoxonRankSum_a_EUD('bs--',2);
        h_um=CGum.fWilcoxonRankSum_a_EUD('ms--',2);
        h_rtog=CGrtog.fWilcoxonRankSum_a_EUD('gs--',2);
        h_comb=CGcomb.fWilcoxonRankSum_a_EUD('ks--',2);
        title('EUD Rank-Sum p-value','FontSize',15);
        set(gca,'box','on');
        xlabel('log_1_0(a)');
        ylabel('p-value');
        h_lgnd=legend([h_msk h_nki h_rtog h_um h_comb],'MSK','NKI','RTOG','UMich','COMB','Location','Best');
        set(h_lgnd,'FontSize',15);
        if do_print, print_fig(gcf,fig_loc,'rank_sum','png');end;
        
        
end
toc;