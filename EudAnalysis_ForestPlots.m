function EudAnalysis_ForestPlots
tic; %close all;

do_print = false;

%% Historical Data

% msk_nki_a = (-1, 0.02]


%fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_meta.mat';

%fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_meta.mat';
%fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_fine_meta.mat';

%% can't load
%fn = 'C:\Documents and
%Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';

fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';

%fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_med_EUD_fine_meta.mat';
%fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat';


load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');
 % log10(n) correction for figures
LymanN = log10(CGmsk.mLymanN);
    
CGmsk.mLymanN = LymanN;
CGnki.mLymanN = LymanN;
CGum.mLymanN = LymanN;
CGrtog.mLymanN = LymanN;
CGcomb.mLymanN = LymanN;
    
msk = CGmsk.mGrp;
nki = CGnki.mGrp;
umich = CGum.mGrp;
rtog = CGrtog.mGrp;
comb = CGcomb.mGrp;

CGs = {CGmsk CGnki CGrtog CGum CGcomb};
CGgrps = {msk nki rtog umich  comb};
groups = {'MSK' 'NKI' 'RTOG' 'UMich' 'COMB'};

logas = cell(length(CGgrps),1); % each cell [loga, [-,+ 68 CL], [-,+ 95 CL]]
weights = zeros(length(CGgrps),1);
comb_loga=0;
comb_68L=0;
comb_68H=0;
comb_95L=0;
comb_95H=0;
comb_wgt=0;
loga_95cl_low=-inf(length(CGgrps),1);
loga_95cl_high=-inf(length(CGgrps),1);
for i=1:length(CGs)
    OCobj = CGs{i};
    % find min loga

    
    
    [mx_loga,loc] = max(OCobj.mLymanGrid.loglikelihood(:));
    [~,~,a_loc] = ind2sub(size(OCobj.mLymanGrid.loglikelihood),loc);
    loga = -OCobj.mLymanN(a_loc);
       
    
    
    %% loga CLs
    llmx = zeros(size(OCobj.mLymanN)); %log likelihood's max for each log10(a)
    for kk = 1:length(llmx)
            ll = OCobj.mLymanGrid.loglikelihood(:,:,kk);
            llmx(kk) = max(ll(:));
    end
    
    llmx = llmx / (OCobj.mNumInGrp-2);
    mx_loga = mx_loga / (OCobj.mNumInGrp-2);
            
    lowlog68 = mx_loga-0.5*1/(OCobj.mNumInGrp-2);
    lowlog95 = mx_loga-0.5*(1.96*2)/(OCobj.mNumInGrp-2);
             
    [~, CI68loga] = ConfidenceInterval(OCobj.mLymanN,flipdim(llmx,1), lowlog68);
    [~, CI95loga] = ConfidenceInterval(OCobj.mLymanN,flipdim(llmx,1), lowlog95);
    CI95loga(CI95loga==-Inf)=-1;
    CI68loga(CI68loga==-Inf)=-1;
    logas{i} = {loga,CI68loga,CI95loga};    

    
    
    %cur_wgt = 1/abs(CI95loga(2)-CI95loga(1));
    %tmp
    cur_wgt = OCobj.mNumInGrp;
    weights(i) = cur_wgt;
    comb_loga = comb_loga+loga*cur_wgt;
    
    if i<length(CGs) % don't add comb
        comb_wgt = comb_wgt+cur_wgt;
    end
    
    comb_68L = comb_68L + (CI68loga(1)*cur_wgt)^2;
    comb_68H = comb_68H + (CI68loga(2)*cur_wgt)^2;
    
    comb_95L = comb_95L + (CI95loga(1)*cur_wgt)^2;
    comb_95H = comb_95H + (CI95loga(2)*cur_wgt)^2;
    
    disp(groups{i});        
    disp(['68% log(a) CIs: ',...
        num2str(loga),...
        ' [',...
        num2str(CI68loga(1)),...
        ', ',...
        num2str(CI68loga(2)),']']);
        
    disp(['95% log(a) CIs: ',...
        num2str(loga),...
        ' [',...
        num2str(CI95loga(1)),...
        ', ',...
        num2str(CI95loga(2)),']']);

    loga_95cl_low(i) = CI95loga(1);
    loga_95cl_high(i) = CI95loga(2);    
    
end
% comb_loga = (comb_loga/length(CGs))/comb_wgt;
% err_loga = std([cellfun(@(x) x{1},logas(1:end-1))]);
% err_loga = (err_loga/sqrt(length(CGs)))*1.96
% 
% weights(end)=comb_wgt;
% logas{end}={comb_loga,[0, 0],[comb_loga-err_loga comb_loga+err_loga]};

weights(end)=comb_wgt;

%%  loga
%h_loga=fForestPlotLogA(CGgrps,logas,weights,groups, 'log_{10}(a)');

h_loga=fForestPlotLogA_YIS(CGgrps,logas,weights,groups,...
    loga_95cl_low,loga_95cl_high, '$\log_{10}(a)$');
end


function fig=fForestPlotLogA(cgs,var,wgt,grps,var_str)
    screen_size=get(0,'ScreenSize');
    %ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
    ss_four2three = [0 0 screen_size(3)/2 screen_size(4)/2];
    figure;clf reset;
    set(gcf,'Position',ss_four2three);
    
    y=1:length(var);
    
    means = flipud([cellfun(@(x) x{1},var)]);
    
    mean_95H = flipud([cellfun(@(x) x{3}(2),var)]);
    mean_95L = flipud([cellfun(@(x) x{3}(1),var)]);
    
    
    fig=herrorbar(means,y',abs(mean_95L-means),abs(mean_95H-means));
    set(fig,'Color','k');
    hold on;
    
    wgt=flipud(wgt);
    tot_wgt=wgt(1);
    symbols = {'kd' 'ks' 'ks' 'ks' 'ks'};
    size_factor = [25 50 50 50 50];
    marker_fc = {'w' 'k' 'k' 'k' 'k'};
    for j=y
        msize = (wgt(j)/tot_wgt)*size_factor(j);
            plot(means(j),y(j),symbols{j},'MarkerSize',msize,...
                'MarkerFaceColor',marker_fc{j});
    end
    %set(gca,'YTick',[]);
    %set(gca,'YColor','w');
    set(gca,'box','off');
    hold on;
    xlim([-3 1]);
    ylim([0 max(y)+1]);
    set(gca,'XTick',-3:0.5:1);
    set(gca,'XTickLabel',{[] [] [] [] '-1' '-0.5' '0' '0.5' '1'});
    set(gca,'YTick',1:length(var));
    set(gca,'YTickLabel',fliplr(grps));
    set(gca,'FontSize',14);

    comb_mean = means(1);
    plot(xlim,[max(y)+0.5 max(y)+0.5],'k-')
    plot([comb_mean comb_mean],[0 max(y)+0.5],'r-.');
    plot([0 0],[0 max(y)+0.5],'k--');
    plot([-2.5 -2.5],[0.001 0.5],'w-') % remove tick marks
    plot([-2.0 -2.0],[0.001 0.5],'w-') % remove tick marks
    plot([-1.5 -1.5],[0.001 0.5],'w-') % remove tick marks
    %plot([0 0], [max(y) max(y)+0.5],'w-');
    
    columns = {['Num of',10,'RP/patients'], var_str,'weight',[var_str ' ($95\%$~CL)']};
    text([-2.75 -2 -1.5 -.25], repmat(max(y)+0.75,1,4),columns,'FontSize',14);
    
    num_pts = [cellfun(@(x) length(x),cgs)]; 
    num_rp = [cellfun(@(x) sum(~[x.mFlgCensor]),cgs)];
    frac_rp_pts = num_rp./num_pts;
    str_rp_pts = cell(length(cgs),1);
    for k=1:length(cgs)
        str_rp_pts{k} =...
            [strcat(num2str(num_rp(k)),'/',num2str(num_pts(k))),10,...
            '(',num2str(frac_rp_pts(k)*100,'%2.1f'),'%)'];
    end
    
    text(repmat(-2.75,5,1),y',flipud(str_rp_pts),'FontSize',14);
    
    str_loga = num2str(means,'%2.2f');
    text(repmat(-2,5,1),y',str_loga,'FontSize',14);
    
    str_wgt = strcat(num2str(100*wgt./wgt(1),'%2.1f'),'%');
    text(repmat(-1.5,5,1),y',str_wgt,'FontSize',14);
end


function fig=fForestPlotLogA_YIS(cgs,var,wgt,grps,...
    loga_ci_low,loga_ci_high,var_str)
    screen_size=get(0,'ScreenSize');
    ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
    %ss_four2three = [0 0 screen_size(3)/2 screen_size(4)/2];
    cur_fig=figure;clf reset;
    set(gcf,'Position',ss_four2three);
    
    
    y=1:length(var);
    
    means = flipud([cellfun(@(x) x{1},var)]);
    
    mean_95H = flipud([cellfun(@(x) x{3}(2),var)]);
    mean_95L = flipud([cellfun(@(x) x{3}(1),var)]);
    
    
    fig=herrorbar(means,y',abs(mean_95L-means),abs(mean_95H-means));
    set(fig,'Color','k');
    hold on;
    
    wgt=flipud(wgt);
    tot_wgt=wgt(1);
    symbols = {'kd' 'ks' 'ks' 'ks' 'ks'};
    size_factor = [25 60 60 60 60];
    marker_fc = {'w' 'k' 'k' 'k' 'k'};
    for j=y
        msize = (wgt(j)/tot_wgt)*size_factor(j);
            plot(means(j),y(j),symbols{j},'MarkerSize',msize,...
                'MarkerFaceColor',marker_fc{j});
    end
    %set(gca,'YTick',[]);
    %set(gca,'YColor','w');
    set(gca,'box','off');
    hold on;
    xlim([-3.5 1]);
    ylim([0 max(y)+1]);
    set(gca,'XTick',-3.5:0.5:1);
    set(gca,'XTickLabel',{[] [] [] [] [] '-1' '-0.5' '0' '0.5' '1'});
    set(gca,'YTick',1:length(var));
    set(gca,'YTickLabel',fliplr(grps));
    set(gca,'FontSize',14);

    comb_mean = means(1);
    plot(xlim,[max(y)+0.5 max(y)+0.5],'k-')
    plot([comb_mean comb_mean],[0 max(y)+0.5],'r-.');
    plot([0 0],[0 max(y)+0.5],'k--');
    
    plot([-3.5 -3.5],[0.001 0.5],'w-') % remove tick marks
    plot([-3 -3],[0.001 0.5],'w-') % remove tick marks
    plot([-2.5 -2.5],[0.001 0.5],'w-') % remove tick marks
    plot([-2.0 -2.0],[0.001 0.5],'w-') % remove tick marks
    plot([-1.5 -1.5],[0.001 0.5],'w-') % remove tick marks
    %plot([0 0], [max(y) max(y)+0.5],'w-');
    
    columns = {['Incidence',10,'of RP'],[' ' var_str],'[$95\%$~CI]',[var_str ' ($95\%$~CI)']};
    text([-2.75 -2.05 -1.55 -.25], repmat(max(y)+0.75,1,4),columns,'FontSize',16,'interpreter','latex');
    
    str_data = {'Comb', 'UMich','RTOG','NKI','MSK'};
    text(repmat(-3.5,5,1),y',str_data,'FontSize',18,'interpreter','latex');
    
    num_pts = [cellfun(@(x) length(x),cgs)]; 
    num_rp = [cellfun(@(x) sum(~[x.mFlgCensor]),cgs)];
    frac_rp_pts = num_rp./num_pts;
    str_rp_pts = cell(length(cgs),1);
    for k=1:length(cgs)
        str_rp_pts{k} =...
            ['$',strcat(num2str(num_rp(k)),'/',num2str(num_pts(k))),'$',10,...
            '$(',num2str(frac_rp_pts(k)*100,'%2.1f'),'\%)$'];
    end
    
    text(repmat(-2.75,5,1),y',flipud(str_rp_pts),'FontSize',16,'interpreter','latex');
    
    str_loga = num2str(means,'%2.2f');
    %text(repmat(-2.5,5,1),y',str_loga,'FontSize',14);
    text(repmat(-2,5,1),y',str_loga,'FontSize',16,'interpreter','latex');
    
    %str_a = str2num(str_loga);
    %str_a=10.^str_a;
    %str_a = strcat(num2str(str_a,'%2.2f'));
    %text(repmat(-1.5,5,1),y',str_a,'FontSize',14);
    
    %loga_ci_low = 10.^loga_ci_low;
    str_ci_low = {num2str(loga_ci_low,'%2.2f')};
    
    %loga_ci_high = 10.^loga_ci_high;
    str_ci_high = num2str(loga_ci_high,'%2.2f');
    
    str_loga_ci = cell(length(loga_ci_low),1);
    for j=1:length(loga_ci_low)
        str_loga_ci{j} = ['[',...
            num2str(loga_ci_low(j),'%2.2f'),'-',...
            num2str(loga_ci_high(j),'%2.2f'),']'];
    end
    text(repmat(-1.7,5,1),y',flipud(str_loga_ci),'FontSize',16,'interpreter','latex');

  
    set(gca,'YTick',[]);
    set(gca,'YColor','w');
    
    q_text = {'Heterogeneity','$\rm{P} = 0.10;~\rm{I}^{2}=0.52$'};
    text([-3.25 -2.2], repmat(0.3,2,1),q_text','FontSize',18,'interpreter','latex');
    
    do_print=true;
    fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';

    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'meta_forest'],'-pdf');
        disp(['Saving ',fig_loc,'meta_forest.pdf...']);
    end;
    
    
    % text(repmat(-1.5,5,1),y',str_a,'FontSize',14);
end