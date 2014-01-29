function EudAnalysis_Summary
tic; %close all;
do_boxplots=true;
do_rtog=false;

screen_size=get(0,'ScreenSize');
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';

%fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_meta.mat';
fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_meta.mat';
load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');



msk = CGmsk.mGrp;
msk_comps = ~[msk(:).mFlgCensor];
msk_eud = [msk(:).mEUD];

nki = CGnki.mGrp;
nki_comps = ~[nki(:).mFlgCensor];
nki_eud = [nki(:).mEUD];

umich = CGum.mGrp;
umich_comps = ~[umich(:).mFlgCensor];
umich_eud = [umich(:).mEUD];

rtog = CGrtog.mGrp;
rtog_comps = ~[rtog(:).mFlgCensor];
rtog_eud = [rtog(:).mEUD];

comb = CGcomb.mGrp;
comb_comps = ~[comb(:).mFlgCensor];
comb_eud = [comb(:).mEUD];

disp('n_pts & n_comps & Overall Incidence')
disp(['MSK & ',num2str(length(msk)),' & ',...
    num2str(length(msk(msk_comps))),' & ',...
    num2str(length(msk(msk_comps))/length(msk),3)]);
disp(['NKI & ',num2str(length(nki)),' & ',...
    num2str(length(nki(nki_comps))),' & ',...
    num2str(length(nki(nki_comps))/length(nki),3)]);
disp(['RTOG & ',num2str(length(rtog)),' & ',...
    num2str(length(rtog(rtog_comps))),' & ',...
    num2str(length(rtog(rtog_comps))/length(rtog),3)]);
disp(['UMich & ',num2str(length(umich)),' & ',...
    num2str(length(umich(umich_comps))),' & ',...
    num2str(length(umich(umich_comps))/length(umich),3)]);



%% EUD Boxplots
loga = [1:-0.1:-1];
if do_boxplots
    for n=1:21
        eud_mean_loc=n;
        
        %msk
        msk_eud_cens = msk_eud(:,~msk_comps);
        msk_eud_comps = msk_eud(:,msk_comps);
        msk_mean_cens = msk_eud_cens(eud_mean_loc,:);
        msk_mean_comps = msk_eud_comps(eud_mean_loc,:);
        msk_group = [ones(size(msk_mean_cens'));...
            ones(size(msk_mean_comps'))+1];
        msk_means = cat(1,msk_mean_cens',msk_mean_comps');
        
        %nki
        nki_eud_cens = nki_eud(:,~nki_comps);
        nki_eud_comps = nki_eud(:,nki_comps);
        nki_mean_cens = nki_eud_cens(eud_mean_loc,:);
        nki_mean_comps = nki_eud_comps(eud_mean_loc,:);
        nki_group = [ones(size(nki_mean_cens'));...
            ones(size(nki_mean_comps'))+1];
        nki_means = cat(1,nki_mean_cens',nki_mean_comps');
        
        %umich
        umich_eud_cens = umich_eud(:,~umich_comps);
        umich_eud_comps = umich_eud(:,umich_comps);
        umich_mean_cens = umich_eud_cens(eud_mean_loc,:);
        umich_mean_comps = umich_eud_comps(eud_mean_loc,:);
        umich_group = [ones(size(umich_mean_cens'));...
            ones(size(umich_mean_comps'))+1];
        umich_means = cat(1,umich_mean_cens',umich_mean_comps');
        
        
        %rtog
        rtog_eud_cens = rtog_eud(:,~rtog_comps);
        rtog_eud_comps = rtog_eud(:,rtog_comps);
        rtog_mean_cens = rtog_eud_cens(eud_mean_loc,:);
        rtog_mean_comps = rtog_eud_comps(eud_mean_loc,:);
        rtog_group = [ones(size(rtog_mean_cens'));...
            ones(size(rtog_mean_comps'))+1];
        rtog_means = cat(1,rtog_mean_cens',rtog_mean_comps');
        
        %comb
        comb_eud_cens = comb_eud(:,~comb_comps);
        comb_eud_comps = comb_eud(:,comb_comps);
        comb_mean_cens = comb_eud_cens(eud_mean_loc,:);
        comb_mean_comps = comb_eud_comps(eud_mean_loc,:);
        comb_group = [ones(size(comb_mean_cens'));...
            ones(size(comb_mean_comps'))+1];
        comb_means = cat(1,comb_mean_cens',comb_mean_comps');
        
        
        
        max_ylim = max([max(msk_means) max(nki_means)...
            max(rtog_means) max(umich_means)])+2;
        
        cur_fig=figure(n);clf reset;
        set(gcf,'Position',ss_four2three);
        
        subplot(2,2,1);
        boxplot(msk_means,msk_group);
        set(gca,'XTick',1:2,'XTickLabel',{'Censored','Complication'},'FontSize',14);
        ylabel(['EUD [Gy]'],'FontSize',14);
        ylim([0 max_ylim]);
        title('MSK','FontSize',14);
        
        subplot(2,2,3);
        boxplot(nki_means,nki_group);
        set(gca,'XTick',1:2,'XTickLabel',{'Censored','Complication'},'FontSize',14);
        ylabel(['EUD [Gy]'],'FontSize',14);
        ylim([0 max_ylim]);
        title('NKI','FontSize',14);
        
        subplot(2,2,2);
        boxplot(rtog_means,rtog_group);
        set(gca,'XTick',1:2,'XTickLabel',{'Censored','Complication'},'FontSize',14);
        ylabel(['EUD [Gy]'],'FontSize',14);
        ylim([0 max_ylim]);
        title('RTOG','FontSize',14);
        
        subplot(2,2,4);
        boxplot(umich_means,umich_group);
        set(gca,'XTick',1:2,'XTickLabel',{'Censored','Complication'},'FontSize',14);
        ylabel(['EUD [Gy]'],'FontSize',14);
        ylim([0 max_ylim]);
        title('UMich','FontSize',14);
        
        if loga(n)<0
            fig_str = ['box_euds_loga_m',strrep(num2str(abs(loga(n))),'.','_')];
        else
            fig_str = ['box_euds_loga_',strrep(num2str(abs(loga(n))),'.','_')];
        end
    
        if do_print,
        %print_fig(gcf,fig_loc,fig_str,'png')
        
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'box_euds_loga-',num2str(n)],'-pdf');
            disp(['Saving ',fig_loc,'box_euds_loga-',num2str(n),'.pdf...']);
        end
    
        cur_fig=figure(2*n);clf reset;
        set(gcf,'Position',ss_four2three);
       
        boxplot(comb_means,comb_group);
        set(gca,'XTick',1:2,'XTickLabel',{'Censored','Complication'},'FontSize',14);
        ylabel(['EUD [Gy]'],'FontSize',14);
        ylim([0 max_ylim]);
        title('MSK + NKI + RTOG + UMich','FontSize',14);
        
        if loga(n)<0
            fig_str = ['comb_box_euds_loga_m',strrep(num2str(abs(loga(n))),'.','_')];
        else
            fig_str = ['comb_box_euds_loga_',strrep(num2str(abs(loga(n))),'.','_')];
        end
        
        if do_print,
            %print_fig(gcf,fig_loc,fig_str,'png');end;
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'comb_box_euds_loga-',num2str(n)],'-pdf');
            disp(['Saving ',fig_loc,'comb_box_euds_loga-',num2str(n),'.pdf...']);
        end
    
    end
end

if do_rtog, % plot max doses for RTOG binned in Tx 
   
    rtog_vols = [rtog(:).mVolCum]'; %(num pts X num bins)
    rtog_lq_doses = [rtog(:).mDoseBins_LQ]';
    rtog_phys_doses = [rtog(:).mDoseBins_org]';
    
    rtog_min_vols = rtog_vols<=0;
    rtog_max_lq_doses = zeros(size(rtog_min_vols,1),1);
    rtog_max_phys_doses = zeros(size(rtog_min_vols,1),1);
    
    for i=1:size(rtog_min_vols,1) % for each patient, find max dose
        rtog_min_vol_idx = find(rtog_min_vols(i,:));
        rtog_min_vol_idx = rtog_min_vol_idx(1);
        if rtog_min_vol_idx==length(rtog_vols)
            rtog_min_vol_idx=length(rtog_vols)-1;
        end
        rtog_max_lq_doses(i)=rtog_lq_doses(i,rtog_min_vol_idx);
        rtog_max_phys_doses(i)=rtog_phys_doses(i,rtog_min_vol_idx);
    end
    
    rtog_txs = [rtog(:).mDoseTx]';
    
    figure;
    set(gcf,'Position',ss_four2three);
    boxplot(rtog_max_lq_doses,rtog_txs);
    set(gca,'FontSize',14);
    set(findobj(gca,'Type','text'),'FontSize',14);
    title('RTOG NTD Dose (alpha/beta=3)','FontSize',16);
    xlabel('Prescription Dose [Gy]','FontSize',15);
    ylabel('Max dose [Gy]','FontSize',15);
    ylim([50 150]);
    
    figure;
    set(gcf,'Position',ss_four2three);
    boxplot(rtog_max_phys_doses,rtog_txs);
    set(gca,'FontSize',14);
    set(findobj(gca,'Type','text'),'FontSize',14);
    title('RTOG Physical Dose (alpha/beta=Inf)','FontSize',16);
    xlabel('Prescription Dose [Gy]','FontSize',15);
    ylabel('Max dose [Gy]','FontSize',15);
    ylim([50 150]);
    %% UM
    
    umich_vols = [umich(:).mVolCum]'; %(num pts X num bins)
    if isempty(umich_vols)
        for j=1:length(umich)
            umich(j) = umich(j).fDiff2Cum();
        end
        umich_vols = [umich.mVolCum]';
    end
    
    
    umich_lq_doses = [umich(:).mDoseBins_LQ]';
    umich_phys_doses = [umich(:).mDoseBins_org]';
    
    umich_min_vols = umich_vols<=0;
    umich_max_lq_doses = zeros(size(umich_min_vols,1),1);
    umich_max_phys_doses = zeros(size(umich_min_vols,1),1);
    
    for i=1:size(umich_min_vols,1) % for each patient, find max dose
        umich_min_vol_idx = find(umich_min_vols(i,:));
        umich_min_vol_idx = umich_min_vol_idx(1);
        if umich_min_vol_idx==length(umich_vols)
            umich_min_vol_idx=length(umich_vols)-1;
        end
            
        umich_max_lq_doses(i)=umich_lq_doses(i,umich_min_vol_idx);
        umich_max_phys_doses(i)=umich_phys_doses(i,umich_min_vol_idx);
    end
    
    umich_txs = [umich(:).mDoseTx]';
    
    % only one patient at 83.84, group with 84
    umich_txs(umich_txs==83.84)=84;
    
    figure;
    set(gcf,'Position',ss_four2three);
    boxplot(umich_max_lq_doses,umich_txs);
    set(gca,'FontSize',14);
    set(findobj(gca,'Type','text'),'FontSize',14);
    title('UMich NTD Dose (alpha/beta=3)','FontSize',16);
    xlabel('Prescription Dose [Gy]','FontSize',15);
    ylabel('Max dose [Gy]','FontSize',15);
    ylim([50 150]);
   
    figure;
    set(gcf,'Position',ss_four2three);
    boxplot(umich_max_phys_doses,umich_txs);
    set(gca,'FontSize',14);
    set(findobj(gca,'Type','text'),'FontSize',14);
    title('UMich Physical Dose (alpha/beta=Inf)','FontSize',16);
    xlabel('Prescription Dose [Gy]','FontSize',15);
    ylabel('Max dose [Gy]','FontSize',15);
    ylim([50 150]);
end

disp([]);
toc;
end