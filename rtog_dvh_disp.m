function rtog_dvh_disp

tic;
screen_size=get(0,'ScreenSize');
fig_loc='Z:/elw/MATLAB/meta_analy/figures/latest/';

ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

PtInfo = classDataFromXls();
CGobjs = classOutcomeAnalysis();
CGobjs.mLymanN = 10.^(-1:0.1:1)';

fn ='Z:/elw/MATLAB/original_data/RTOG9311/2013_01_25/TEST_RTOG9311_LUNG-GTV_HETERO_DVH.xlsx';

xlsWholeLungInfo=xlsFileRead(fn);

PtInfo.mXlsRaw = xlsWholeLungInfo(4).xlsRaw;
PtInfo.mLabel = 'Sum of differences';
PtInfo = PtInfo.fExtractColData();
flgdiff = PtInfo.mFlg;
diffs = [cell2mat(PtInfo.mData(flgdiff))];

PtInfo.mXlsRaw = xlsWholeLungInfo(4).xlsRaw;
PtInfo.mLabel = 'Lung-GTV Area';
PtInfo = PtInfo.fExtractColData();
flgarea = PtInfo.mFlg;
areas = [cell2mat(PtInfo.mData(flgarea))];

diffs = diffs./areas;

[n1, xout1]=hist(diffs,[min(diffs):(max(diffs)-min(diffs))/100:max(diffs)]);


f(1)=figure(1);clf reset;
set(gcf,'Position',ss_four2three);
plot(xout1,n1,'LineWidth',2);
set(gca,'fontsize',12);
xlim([min(diffs) max(diffs)]);
xlabel('$\sum_{\rm{DVH~bins}}~$((Lung-GTV) - (Lung-PTV))/(Lung-GTV)','Interpreter','LaTex','fontsize',16);
ylabel('Frequency','fontsize',16);
grid on;


[n2, xout2]=hist(diffs,[-0.1:0.01:0.1]);

f(2)=figure(2);clf reset;
set(gcf,'Position',ss_four2three);
plot(xout2(1:end-1),n2(1:end-1),'LineWidth',2);% leave off overflow bin
set(gca,'fontsize',12);
xlabel('$\sum_{\rm{DVH~bins}}~$((Lung-GTV) - (Lung-PTV))/(Lung-GTV)','Interpreter','LaTex','fontsize',16);
ylabel('Frequency','fontsize',16);
grid on;

[n3, xout3]=hist(diffs,[-0.01:0.001:0.05]);

f(3)=figure(3);clf reset;
set(gcf,'Position',ss_four2three);
plot(xout3(2:end-1),n3(2:end-1),'LineWidth',2);
set(gca,'fontsize',12);
xlabel('$\sum_{\rm{DVH~bins}}$~((Lung-GTV) - (Lung-PTV))/(Lung-GTV)','Interpreter','LaTex','fontsize',16);
ylabel('Frequency','fontsize',16);
grid on;

return;

% display dvhs
PtInfoGTV = classDataFromXls();
PtInfoGTV.mXlsRaw = xlsWholeLungInfo(2).xlsRaw;
PtInfoGTV.mLabel = 'DOSE';
PtInfoGTV = PtInfoGTV.fExtractRowData();
flgdose = PtInfoGTV.mFlg;
dosebins = [cell2mat(PtInfoGTV.mData(flgdose));0];
mNumDoseBins = length(dosebins);


PtInfoGTV.mLabel = 'DOSE'; % find the pt Ids
PtInfoGTV = PtInfoGTV.fExtractColData();
flgPtIdGTV = PtInfoGTV.mFlg; % pt with ID number are flaged
mFullPtIdsGTV = PtInfoGTV.mData;

PtInfoPTV = classDataFromXls();
PtInfoPTV.mXlsRaw = xlsWholeLungInfo(1).xlsRaw;
PtInfoPTV.mLabel = 'DOSE'; % find the pt Ids
PtInfoPTV = PtInfoPTV.fExtractColData();
flgPtIdPTV = PtInfoPTV.mFlg; % pt with ID number are flaged
mFullPtIdsPTV = PtInfoPTV.mData;



mFullPtIds = mFullPtIdsGTV(flgPtIdGTV&flgPtIdPTV);


mPtIds = unique(mFullPtIds);

disp(['To be excluded?');
mPtIds(diffs<=0.01)

[sorted_diffs,idx_diffs] = sort(diffs);
sorted_mPtIds = mPtIds(idx_diffs);

%mPtIds = mPtIds(diffs>1000 & diffs<5000);
mPtIds = sorted_mPtIds;
diffs = sorted_diffs;

mNumPts = length(mPtIds);

%curDiffs = diffs(diffs>1000 & diffs<5000);
dosebins=dosebins(1:end-1);

for i=1:mNumPts
    cur_pt = mPtIds(i);
    PtInfoGTV.mLabel = cur_pt;
    PtInfoGTV = PtInfoGTV.fExtractRowData();
    flgdvhgtv = PtInfoGTV.mFlg;
    gtv_cvol = cell2mat(PtInfoGTV.mData(flgdvhgtv));
    gtv_cvol=gtv_cvol./gtv_cvol(1);
    gtv_cvol(gtv_cvol<=0)=NaN;
    
    PtInfoPTV.mLabel = cur_pt;
    PtInfoPTV = PtInfoPTV.fExtractRowData();
    flgdvhptv = PtInfoPTV.mFlg;
    ptv_cvol = cell2mat(PtInfoPTV.mData(flgdvhptv));
    ptv_cvol=ptv_cvol./ptv_cvol(1);
    ptv_cvol(ptv_cvol<=0)=NaN;
    
    cur_fig=figure(100+i);clf reset;
    set(gcf,'Position',ss_four2three);
    hold on;
    h(1)=plot(dosebins,ptv_cvol,'LineWidth',3);
    h(2)=plot(dosebins,gtv_cvol,'r','LineWidth',3);
    hold off;
    grid on;
    xlabel('Dose [Gy]','FontSize',14);
    ylabel('Volume [cc]','FontSize',14);
    lgnd=legend(h,'Lung-GTV','Lung-PTV');
    set(lgnd,'FontSize',18);
    delta_v = '\Delta ';
    delta_v= [delta_v, sprintf('V = %2.2d cc',diffs(i))];
    title([strrep(cur_pt,'_','\_'),' ',delta_v],'FontSize',18);
    
   fig_str=sprintf('DVH_%d',i);
   print_fig(gcf,fig_loc,fig_str,'png')
   %print(gcf,print_loc,fig_name,'png');
end

end