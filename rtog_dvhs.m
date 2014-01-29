function rtog_dvhs
fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_EUD_meta.mat';
%fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_Homo_EUD_meta_comp_lt6m.mat';
load(fn_rtog,'CGobjs'); CGrtog = CGobjs(1);


%# DVHS
f = [CGrtog.mGrp.mFlgCensor];

avg_dosebins = [CGrtog.mGrp(1).mDoseBins_org];
avg_dosebins = avg_dosebins(1:end-1);
avg_cvol = [CGrtog.mGrp.mVolCum]';
avg_cens_cvol = mean(avg_cvol(f,:));
avg_cens_cvol = avg_cens_cvol(1:end-1);
std_cens_cvol = std(avg_cvol(f,:));
std_cens_cvol = std_cens_cvol(1:end-1);

avg_comp_cvol = mean(avg_cvol(~f,:));
avg_comp_cvol = avg_comp_cvol(1:end-1);
std_comp_cvol = std(avg_cvol(~f,:));
std_comp_cvol = std_comp_cvol(1:end-1);



figure(1); clf reset; hold on; grid on;
hold on;
h(1)=plot(avg_dosebins,avg_comp_cvol,'r','LineWidth',2);
h(2)=plot(avg_dosebins,avg_comp_cvol+std_comp_cvol,'--r','LineWidth',0.5);
plot(avg_dosebins,avg_comp_cvol-std_comp_cvol,'--r','LineWidth',0.5);

h(3)=plot(avg_dosebins,avg_cens_cvol,'b','LineWidth',2);
h(4)=plot(avg_dosebins,avg_cens_cvol+std_cens_cvol,'--b','LineWidth',0.5);
plot(avg_dosebins,avg_cens_cvol-std_cens_cvol,'--b','LineWidth',0.5);

legend(h,'Avg. no comp','Central 68%','Avg. with comp','Central 68%');
ylim([0,1]);
xlabel('Dose [Gy]','FontSize',15);
ylabel('Volume Fraction','FontSize',15);
hold off;

figure(2); clf reset; hold on; % grid on;
% DVHs of censored patients
for k=1:length(f),
    if f(k), %# cens
        dvh_color = 'b';
    else
        dvh_color = 'r';
    end
    
    dosebins = CGrtog.mGrp(k).mDoseBins_org;
    volcum = CGrtog.mGrp(k).mVolCum;
    plot(dosebins(1:end-1), volcum(1:end-1),dvh_color);
    
end
set(gca,'xminortick','on','yminortick','on');
xlabel('Dose [Gy]','FontSize',15);
ylabel('Volume Fraction','FontSize',15);





mLymanN = -1:0.1:1;

avg_eud = [CGrtog.mGrp.mEUD]';

avg_cens_eud = mean(avg_eud(f,:));
std_cens_eud = std(avg_eud(f,:));

avg_comp_eud = mean(avg_eud(~f,:));
std_comp_eud = std(avg_eud(~f,:));


figure(3); clf reset; hold on; grid on;
hold on;
% h(1)=plot(-mLymanN,avg_comp_eud,'r','LineWidth',2);
% h(2)=plot(-mLymanN,avg_comp_eud+std_comp_eud,'--r','LineWidth',0.5);
% plot(-mLymanN,avg_comp_eud-std_comp_eud,'--r','LineWidth',0.5);
% 
% h(3)=plot(-mLymanN,avg_cens_eud,'b','LineWidth',2);
% h(4)=plot(-mLymanN,avg_cens_eud+std_cens_eud,'--b','LineWidth',0.5);
% plot(-mLymanN,avg_cens_eud-std_cens_eud,'--b','LineWidth',0.5);
h(1)=plot(avg_comp_eud,mLymanN,'-r*','LineWidth',2);
h(2)=plot(avg_comp_eud+std_comp_eud,mLymanN,'--r','LineWidth',0.5);
plot(avg_comp_eud-std_comp_eud,mLymanN,'--r','LineWidth',0.5);

h(3)=plot(avg_cens_eud,mLymanN,'-b*','LineWidth',2);
h(4)=plot(avg_cens_eud+std_cens_eud,mLymanN,'--b','LineWidth',0.5);
plot(avg_cens_eud-std_cens_eud,mLymanN,'--b','LineWidth',0.5);

set(gca,'YTickLabel',1:-0.2:-1);
legend(h,'Avg. no comp','Central 68%','Avg. with comp','Central 68%',...
    'Location','NorthEast');
ylabel('log_1_0(a)','FontSize',15);
xlabel('gEUD','FontSize',15);


hold off;



figure(4); clf reset; hold on; % grid on;
% DVHs of censored patients
for k=1:length(f),
    if f(k), %# cens
        dvh_color = 'b';
    else
        dvh_color = 'r';
    end
    
    euds = CGrtog.mGrp(k).mEUD';
    plot(euds,mLymanN,dvh_color);
    
end
set(gca,'YTickLabel',1:-0.2:-1);
set(gca,'xminortick','on','yminortick','on');
ylabel('log_1_0(a)','FontSize',15);
xlabel('gEUD [Gy]','FontSize',15);



% figure(3); clf reset; hold on; % grid on;
% f = [CGrtog.mGrp.mFlgCensor];
% % DVHs of censored patients
% for k=1:length(f),
%     if f(k), %# cens
%         dvh_color = 'b';
%     else
%         dvh_color = 'r';
%     end
%     dosebins = CGrtog.mGrp(k).mDoseBins_org;
%     voldiff = CGrtog.mGrp(k).mVolDiff;
%     plot(dosebins(1:end-1), voldiff(1:end-1),dvh_color);
% end
% ylim([0,0.15]);
% set(gca,'xminortick','on','yminortick','on');
% xlabel('Dose [Gy]','FontSize',15);
% ylabel('Volume Fraction','FontSize',15);

end