function EudAnalysis_Correlations
tic; %close all;

screen_size=get(0,'ScreenSize');
ss_two2two = [screen_size(3)/2 0 screen_size(4) screen_size(4)];
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;

fig_loc = 'Z:/elw/MATLAB/meta_analy/figures/latest/';
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


%% EUD correlations
loga = [1:-0.1:-1];
msk_eud_corrs = corr(msk_eud',msk_eud').^2;
nki_eud_corrs = corr(nki_eud',nki_eud').^2;
rtog_eud_corrs = corr(rtog_eud',rtog_eud'.^2);
umich_eud_corrs = corr(umich_eud',umich_eud').^2;

figure;clf reset;
set(gcf,'Position',ss_four2three);

subplot(2,2,1);
imagesc(loga,loga,msk_eud_corrs);
set(gca,'YDir','normal');
colorbar;
caxis([0 1]);
set(gca,'FontSize',12);
ylabel(['log_{10}(a)'],'FontSize',14);
xlabel(['log_{10}(a)'],'FontSize',14);
title('MSK r^2','FontSize',14);

subplot(2,2,2);
imagesc(loga,loga,rtog_eud_corrs);
set(gca,'YDir','normal');
colorbar;
caxis([0 1]);
set(gca,'FontSize',12);
ylabel(['log_{10}(a)'],'FontSize',14);
xlabel(['log_{10}(a)'],'FontSize',14);
title('RTOG r^2','FontSize',14);

subplot(2,2,3);
imagesc(loga,loga,umich_eud_corrs);
set(gca,'YDir','normal');
colorbar;
caxis([0 1]);
set(gca,'FontSize',12);
ylabel(['log_{10}(a)'],'FontSize',14);
xlabel(['log_{10}(a)'],'FontSize',14);
title('UMich r^2','FontSize',14);

subplot(2,2,4);
imagesc(loga,loga,nki_eud_corrs);
set(gca,'YDir','normal');
colorbar;
caxis([0 1]);
set(gca,'FontSize',12);
ylabel(['log_{10}(a)'],'FontSize',14);
xlabel(['log_{10}(a)'],'FontSize',14);
title('NKI r^2','FontSize',14);


end