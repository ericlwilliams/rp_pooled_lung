function EudAnalysis_DVHs
tic; %close all;
screen_size=get(0,'ScreenSize');

ss_four2three = [0 0 screen_size(3)/2 screen_size(4)/2];

do_print = true;

fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';
%fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_meta.mat';
fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_meta.mat';
load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');

%fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_EUD_meta.mat';
%fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_Homo_EUD_meta_comp_lt6m.mat';
%load(fn_rtog,'CGobjs'); CGrtog = CGobjs(1);
CGs = {CGrtog CGum CGnki};
titles = {'RTOG' 'UMich' 'NKI'};
for i=1:length(CGs)
    
    CGcur = CGs{i};
    %# DVHS
    f = [CGcur.mGrp.mFlgCensor];
    
    avg_dosebins = [CGcur.mGrp(1).mDoseBins_LQ];
    avg_dosebins = avg_dosebins(1:end-1);
    avg_cvol = [CGcur.mGrp.mVolCum]';
    if isempty(avg_cvol)
        for j=1:length(CGcur.mGrp)
            CGcur.mGrp(j) = CGcur.mGrp(j).fDiff2Cum();
        end
        avg_cvol = [CGcur.mGrp.mVolCum]';
    end
    
    avg_cens_cvol = mean(avg_cvol(f,:));
    avg_cens_cvol = avg_cens_cvol(1:end-1);
    std_cens_cvol = std(avg_cvol(f,:));
    std_cens_cvol = std_cens_cvol(1:end-1);
    
    avg_comp_cvol = mean(avg_cvol(~f,:));
    avg_comp_cvol = avg_comp_cvol(1:end-1);
    std_comp_cvol = std(avg_cvol(~f,:));
    std_comp_cvol = std_comp_cvol(1:end-1);
    
    
    
    figure(1); hold on;
    if i==1, set(gcf,'Position',ss_four2three);end;
    
    
    subplot(2,2,i+1);
    hold on;
    h(1)=plot(avg_dosebins,avg_comp_cvol,'r','LineWidth',2);
    h(2)=plot(avg_dosebins,avg_comp_cvol+std_comp_cvol,'--r','LineWidth',0.5);
    plot(avg_dosebins,avg_comp_cvol-std_comp_cvol,'--r','LineWidth',0.5);
    
    h(3)=plot(avg_dosebins,avg_cens_cvol,'b','LineWidth',2);
    h(4)=plot(avg_dosebins,avg_cens_cvol+std_cens_cvol,'--b','LineWidth',0.5);
    plot(avg_dosebins,avg_cens_cvol-std_cens_cvol,'--b','LineWidth',0.5);
    
    legend(h,'Avg. with comp','Central 68%','Avg. no comp','Central 68%');
    ylim([0 1]);
    xlim([0 125]);
    xlabel('Dose [Gy]','FontSize',15);
    ylabel('Volume Fraction','FontSize',15);
    title(titles{i},'FontSize',15);
    grid on;
    

    
    figure(2);hold on; % grid on;
    if i==1
        set(gcf,'Position',ss_four2three);
    end
    subplot(2,2,i+1);
    % DVHs of censored patients
    hold on;
    for k=1:length(f),
        if f(k), %# cens
            dvh_color = 'b';
        else
            dvh_color = 'r';
        end
        
        dosebins = CGcur.mGrp(k).mDoseBins_LQ;
        volcum = CGcur.mGrp(k).mVolCum;
        plot(dosebins(1:end-1), volcum(1:end-1),dvh_color);
        
    end
    xlim([0 125]);
    ylim([0 1]);
    title(titles{i},'FontSize',15);
    set(gca,'xminortick','on','yminortick','on');
    xlabel('Dose [Gy]','FontSize',15);
    ylabel('Volume Fraction','FontSize',15);
    
    
    
    
   
    
    
end

CGs = {CGmsk CGrtog CGum CGnki};
titles = {'MSK' 'RTOG' 'UMich' 'NKI'};
names = {'msk' 'rtog' 'um' 'nki'};


for k=1:length(CGs)
    
    mLymanN = -1:0.1:1;
    
    CGcur = CGs{k};
    %# DVHS
    f = [CGcur.mGrp.mFlgCensor];
    
    avg_eud = [CGcur.mGrp.mEUD]';
    
    avg_cens_eud = mean(avg_eud(f,:));
    std_cens_eud = std(avg_eud(f,:));
    
    avg_comp_eud = mean(avg_eud(~f,:));
    std_comp_eud = std(avg_eud(~f,:));
    
    
    fig_avg_euds=figure(k*10); hold on; grid on;
    set(gcf,'Position',ss_four2three);
    %if k==1, set(gcf,'Position',ss_four2three);end;
    %subplot(2,2,k);
    hold on;
    
    h(1)=plot(avg_comp_eud,mLymanN,'-r*','LineWidth',2);
    h(2)=plot(avg_comp_eud+std_comp_eud,mLymanN,'--r','LineWidth',0.5);
    plot(avg_comp_eud-std_comp_eud,mLymanN,'--r','LineWidth',0.5);
    
    h(3)=plot(avg_cens_eud,mLymanN,'-b*','LineWidth',2);
    h(4)=plot(avg_cens_eud+std_cens_eud,mLymanN,'--b','LineWidth',0.5);
    plot(avg_cens_eud-std_cens_eud,mLymanN,'--b','LineWidth',0.5);
    
    ylim([-1 1]);
    xlim([0 85]);
    
    title('Average gEUDs','FontSize',18);
    set(gca,'YTickLabel',1:-0.5:-1);
    legend(h,'Avg. with comp','Central 68%','Avg. no comp','Central 68%',...
        'Location','NorthEast');
    ylabel('log_1_0(a)','FontSize',18);
    xlabel('gEUD [Gy]','FontSize',18);
    grid on;
  
      
    if do_print,
        set(fig_avg_euds,'Color','w');
        export_fig(fig_avg_euds,[fig_loc,names{k},'_avg_euds'],'-pdf');
        disp(['Saving ',fig_loc,names{k},'_avg_euds.pdf...']);
    end;
    
    % spaghetti euds
    
    mLymanN = -1:0.1:1;
    
    figure(4); hold on; % grid on;
    if k==1, set(gcf,'Position',ss_four2three);end;
    subplot(2,2,k);
    hold on;
    % DVHs of censored patients
    first_comp=true;
    first_cens=true;
    for l=1:length(f),
        if f(l), %# cens
            dvh_color = 'b';
        else
            dvh_color = 'r';
        end
        
        euds = CGcur.mGrp(l).mEUD';
        if first_comp && ~f(l)
            g(1)=plot(euds,mLymanN,dvh_color);
            first_comp=false;
        elseif first_cens && f(l)
            g(2)=plot(euds,mLymanN,dvh_color);
            first_cense=false;
        else
            plot(euds,mLymanN,dvh_color);
        end
                      
    end
    ylim([-1 1]);
    xlim([0 85]);
    title(titles{k},'FontSize',15);
    set(gca,'YTickLabel',1:-0.5:-1);
    set(gca,'xminortick','on','yminortick','on');
    legend(g,'Complication','No Complication',...
        'Location','NorthEast');
    
    ylabel('log_1_0(a)','FontSize',15);
    xlabel('gEUD [Gy]','FontSize',15);

    
end
end