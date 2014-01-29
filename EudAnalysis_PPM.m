function EudAnalysis_PPM
% PPM - perfusion parallel model
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];
  
do_print = true;
do_spaghetti = true; % F-speedup

    fig_loc = 'Z:/elw/MATLAB/meta_analy/slides/figures/latest/';
    
tic; %close all;
fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta_ppm.mat';
grp_labels= {'msk','nki','rtog','umich','comb'};

load(fn,'CGrtog','CGnki','CGmsk','CGum','CGcomb');


groups = {CGmsk.mGrp CGnki.mGrp CGrtog.mGrp CGum.mGrp CGcomb.mGrp};
clear CGmsk;
clear CGnki;
clear CGrtog;
clear CGum;
clear CGcomb;

fit_m = 0.426; % 42.6 (40.7-44.6)%
fit_d50 = 28.7; % 28.7 (26.3-31.1) Gy
fit_k = 2.2; % 2.2 (1.8-2.5)

fit_hi_m = 0.446;
fit_hi_d50 = 26.3;
fit_hi_k = 1.8;

fit_lo_m = 0.407;
fit_lo_d50 = 31.1;
fit_lo_k = 2.5;


% NTD with a2b = 2 already calculated

dvhs = cell(length(groups),1);
dprs = cell(length(groups),1);
perfs = cell(length(groups),1);
perfs_hi =  cell(length(groups),1);
perfs_lo =  cell(length(groups),1);

pttotals = cell(length(groups),1);
ptcomps = cell(length(groups),1);

probit_llhds = inf(length(groups),1);
probit_llhds_pdf = inf(length(groups),1);
probit_pvals = inf(length(groups),1);


probit_stats = cell(length(groups),1);
probit_bs = cell(length(groups),1);
probit_hi_bs = cell(length(groups),1);
probit_lo_bs = cell(length(groups),1);

max_perf=0;
max_ntd=0;



for i=1:length(groups)

    if do_spaghetti, 
        f_dvhs=figure(99+i); clf reset; 
        set(f_dvhs,'Position',ss_four2three);
    end

    dvhs{i} = {groups{i}.mDoseBins_LQ};
        
    cur_dprs = cell(length(dvhs{i}),1);
    cur_hi_dprs = cell(length(dvhs{i}),1);
    cur_lo_dprs = cell(length(dvhs{i}),1);

    
    cur_perfs = inf(length(dvhs{i}),1);
    cur_hi_perfs = inf(length(dvhs{i}),1);
    cur_lo_perfs = inf(length(dvhs{i}),1);
    
    for j=1:length(dvhs{i})
        cur_dvh = dvhs{i}{j};
        % convert to perfusion-reduction-histograms
        cur_dprs{j} = fit_m./(1+(fit_d50./cur_dvh).^(fit_k));
        
        %below is confusing but correct
        %higher values of perfusion loss, lead to lower response function
        cur_hi_dprs{j} = fit_lo_m./(1+(fit_lo_d50./cur_dvh).^(fit_lo_k));
        cur_lo_dprs{j} = fit_hi_m./(1+(fit_hi_d50./cur_dvh).^(fit_hi_k));
        
        CGind = groups{i}(j);
        is_comp = ~[groups{i}.mFlgCensor];
        is_comp = is_comp(j);
        %CGind=CGind.fDiff2Cum();
        vols = CGind.mVolDiff;
        vols = vols./sum(vols); %normalize
        
        % calculate pt perfusion
        cur_perfs(j) = sum(cur_dprs{j}.*vols);     
        cur_hi_perfs(j) = sum(cur_hi_dprs{j}.*vols);     
        cur_lo_perfs(j) = sum(cur_lo_dprs{j}.*vols);     

        if do_spaghetti,
            figure(99+i);hold on;
            %rebin_x = [1:150];
            %[new_x,new_y]=rebin(cur_dvh,vols.*100,rebin_x);
            %new_y=(new_y./100);
        
            if is_comp
                plot(cur_dvh,vols,'r','LineWidth',1.5);
            else
                plot(cur_dvh,vols,'b');
            end
        end
        
    end
    
    if max(cur_perfs)>max_perf
        max_perf=max(cur_perfs);
    end
    
    
       
    dprs{i} = cur_dprs;
    perfs{i} = cur_perfs;
    perfs_hi{i} = cur_hi_perfs;
    perfs_lo{i} = cur_lo_perfs;
    
    %% Do some fits
    pttotals{i} = ones(size(groups{i},1),1);    
    ptcomps{i} = ones(size(groups{i},1),1);     
    ptcomps{i}([groups{i}.mFlgCensor])=0;   

    [probit_bs{i},~,probit_stats{i}]=glmfit(perfs{i},[ptcomps{i} pttotals{i}],'binomial','link','probit');  
    probit_pvals(i) = probit_stats{i}.p(2);
    
    %CIs
    probit_hi_bs{i}=glmfit(perfs_hi{i},[ptcomps{i} pttotals{i}],'binomial','link','probit');  

    probit_lo_bs{i}=glmfit(perfs_lo{i},[ptcomps{i} pttotals{i}],'binomial','link','probit');  

    cur_pr=normcdf((perfs{i}.*probit_bs{i}(2))+probit_bs{i}(1));
    cur_pr(~ptcomps{i}) = 1-cur_pr(~ptcomps{i}); % non-complication patients
    cur_pr = log(cur_pr); % log likelihood of each patients
    cur_ll = sum(cur_pr); % loglikelihood of all

    probit_llhds(i) = cur_ll;
    probit_llhds_pdf(i) = cur_ll/(length(ptcomps{i})-length(probit_stats{i}.p));
    
    if do_spaghetti,
        cur_fig=figure(99+i);
        ylim([0 0.25]);
        xlim([0 120]);
        set(gca,'FontSize',18);
        xlabel('NTD (Gy_{3})','FontSize',20);
        ylabel('Volume (%)','FontSize',20);
        
        if do_print,
            set(cur_fig,'Color','w');
            export_fig(cur_fig,[fig_loc,'h_',grp_labels{i},'_dvhs'],'-pdf');
            disp(['Saving ',fig_loc,'h_',grp_labels{i},'_dvhs.pdf...']);
         end

    end
end



disp(['      ',grp_labels]);
disp(['LLHD: ',num2str(probit_llhds')]);
disp(['LLHD/df: ',num2str(probit_llhds_pdf')]);
disp(['AIC: ',num2str(2*2-2.*probit_llhds')]);

%% log-likelihood ratio test
geud_stats = -2.*[-26.23 -36.20 -24.97 -22.08 -120.39];%from EudAnalysis_Meta_disp.m
geud_mat = repmat(geud_stats,5,1);
perf_stats = -2.*probit_llhds;
perf_mat = repmat(perf_stats,1,5);

%% diff_mat
% | (perf_msk-geud_msk) (perf_msk-geud_nki) (perf_msk-geud_rtog) ...|
% | (perf_nki-geud_msk) (perf_nki-geud_nki) (perf_nki-geud_rtog) ...|
% | (perf_rtog-geud_msk) (perf_rtog-geud_nk) (perf_rtog-geud_rtog) ...|
diff_mat = perf_mat-geud_mat;

llr_pvals = 1-chi2cdf(diff_mat,1);

%ignoring off-diagnols
disp(['LLR p: ',num2str(diag(llr_pvals)')]);


%  uni_llr_stats = diag(mv_llr_stats);
%     
%     llr_pvals1 = zeros(size(mv_llr_stats));
%     llr_pvals2 = zeros(size(mv_llr_stats));
%     for i=1:8
%         cur_llr_stats = uni_llr_stats(i)-mv_llr_stats(i,:);
%         llr_pvals1(i,:) = 1-chi2cdf(cur_llr_stats,1); % 1 additional variable
%         cur_llr_stats2 = uni_llr_stats(i)-mv_llr_stats(:,i);
%         llr_pvals2(:,i) = 1-chi2cdf(cur_llr_stats2,1); % 1 additional variable
%     end
%     llr_pvals = max(llr_pvals1,llr_pvals2);
%     llr_pvals(llr_pvals==1)=NaN;
%     
    
    


%% Perfusion model

ntd_x = 1:5:250;
perf_y = fit_m./(1+(fit_d50./ntd_x).^(fit_k)); 
perf_hi_y = fit_hi_m./(1+(fit_hi_d50./ntd_x).^(fit_hi_k)); 
perf_lo_y = fit_lo_m./(1+(fit_lo_d50./ntd_x).^(fit_lo_k)); 

cur_fig=figure(1); clf reset; 
set(cur_fig,'Position',ss_four2three);
h_perf=plot(ntd_x,perf_y,'k-','LineWidth',2);hold on;
h_hi=plot(ntd_x,perf_hi_y,'k-.','LineWidth',1);
h_lo=plot(ntd_x,perf_lo_y,'k-.','LineWidth',1);
set(gca,'FontSize',18);
xlabel('NTD (Gy_3)','FontSize',20);
ylabel('Perfusion Reduction (%)','FontSize',20);
xlim([0 250]);
ylim([0 0.45]);
lgnd1=legend([h_perf h_hi],'Logistic model','95% CI','Location','Best');
set(lgnd1,'FontSize',18);

ax1=gca;
set(ax1,'xlim',[0 250],'ylim',[0 0.45],...
    'XColor','k','YColor','k');

if do_print,
    set(cur_fig,'Color','w');
    export_fig(cur_fig,[fig_loc,'perf_red_fit'],'-pdf');
   disp(['Saving ',fig_loc,'perf_red_fit.pdf...']);
end

for i=1:length(groups)
    
       % perfusion ranksum
    cur_fig=figure(i+1); clf reset; 
    set(cur_fig,'Position',ss_four2three);
    
    [sort_perfs,sort_perf_idx] = sort(perfs{i});
    sort_ptcomps = logical(ptcomps{i}(sort_perf_idx));
    perf_x = 1:length(sort_perfs);
    h_cens=plot(perf_x(~sort_ptcomps),sort_perfs(~sort_ptcomps),'bo','MarkerSize',20);hold on;
    h_comp=plot(perf_x(sort_ptcomps),sort_perfs(sort_ptcomps),'ro','MarkerSize',20,'MarkerFaceColor','r');
    
    ylim([0 max_perf]);
    xlim([0 max(perf_x)+1]);
    lgnd=legend([h_cens h_comp],'No complication','Complication','Location','SouthEast');
    set(lgnd,'FontSize',18);
    set(gca,'FontSize',18);
    xlabel(['Rank'],'FontSize',20);
    ylabel(['Perfusion Reduction (%)'],'FontSize',20);
    
    [rs_p,~,rs_st]=ranksum(sort_perfs(sort_ptcomps),sort_perfs(~sort_ptcomps));
    
    text(0.1,0.8,...
        ['Ranksum:',10,...
        '~~~z = ',num2str(rs_st.zval,4),10,...
        '~~~p = ',num2str(rs_p,2)],...
        'FontSize',25,'Units','normalized','interpreter','latex')
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'ranked_',grp_labels{i},'_perfs'],'-pdf');
        disp(['Saving ',fig_loc,'ranked_',grp_labels{i},'_perfs.pdf...']);
    end
    
     %% profit fits
    cur_fig=figure(i+1+length(groups)); clf reset; 
    set(cur_fig,'Position',ss_four2three);
    
    %cur_perf_x = (min(perfs{:})-0.001):0.001:(max(perfs{:})+0.001);
    cur_perf_x = (0:0.001:(max_perf+0.001));
 
    %probit_y = glmval(probit_bs{i},cur_perf_x,'probit','size',ones(length(cur_perf_x),1));
    %probit_hi_y = glmval(probit_hi_bs{i},cur_perf_x,'probit','size',ones(length(cur_perf_x),1));
    %probit_lo_y = glmval(probit_lo_bs{i},cur_perf_x,'probit','size',ones(length(cur_perf_x),1));
    
    %% testing adding errors in quad
    
    [probit_y,prob_err_lo_y,prob_err_hi_y] =...
        glmval(probit_bs{i},cur_perf_x,'probit',probit_stats{i});

    perf_err_hi_y = glmval(probit_hi_bs{i},cur_perf_x,'probit','size',ones(length(cur_perf_x),1));
    perf_err_lo_y = glmval(probit_lo_bs{i},cur_perf_x,'probit','size',ones(length(cur_perf_x),1));

    perf_err_hi_y = perf_err_hi_y-probit_y;
    perf_err_lo_y = probit_y-perf_err_lo_y;
    
    
    comb_err_hi_y = (((perf_err_hi_y).^2+(prob_err_hi_y).^2).^.5);
    comb_err_lo_y = (((perf_err_lo_y).^2+(prob_err_lo_y).^2).^.5);
            
    plot(cur_perf_x,probit_y,'k.-','LineWidth',3);hold on;
    
    %errors from perfusion reduction fit
    h_perf_err=plot(cur_perf_x,probit_y+perf_err_hi_y,'g--','LineWidth',2);
    plot(cur_perf_x,probit_y-perf_err_lo_y,'g--','LineWidth',2);
    
    %errors from probit fit
    h_prob_err=plot(cur_perf_x,probit_y+prob_err_hi_y,'b--','LineWidth',2);
    plot(cur_perf_x,probit_y-prob_err_lo_y,'b--','LineWidth',2);
    
    % combined errors
    h_comb_err=plot(cur_perf_x,probit_y+comb_err_hi_y,'r--','LineWidth',2);
    plot(cur_perf_x,probit_y-comb_err_lo_y,'r--','LineWidth',2);
    
    
    xlim([min(cur_perf_x) max(cur_perf_x)]);
    ylim([0 0.7]);

    lgnd=legend([h_perf_err h_prob_err h_comb_err],...
            '$\sigma_{\rm{perf}}$',...            
            '$\sigma_{\rm{probit}}$',...
            '$\sqrt{\sigma_{\rm{probit}}^{2}+\sigma_{\rm{perf}}^{2}}$',...
            'location','best');
    set(lgnd,'interpreter','latex');
    set(lgnd,'fontsize',24);
    
    quarts = quantile(perfs{i},3);
    first_quart = logical(perfs{i}<quarts(1));
    second_quart = logical(perfs{i}<quarts(2));
    third_quart = ~second_quart;
    second_quart = logical(second_quart.*~first_quart);
    fourth_quart=logical(perfs{i}>=quarts(3));
    third_quart = logical(third_quart.*~fourth_quart);
    
    meds = [median(perfs{i}(first_quart)),...
            median(perfs{i}(second_quart)),...
            median(perfs{i}(third_quart)),...
            median(perfs{i}(fourth_quart))];
    num_quarts = [sum(first_quart),...
                  sum(second_quart),...
                  sum(third_quart),...  
                  sum(fourth_quart)];
    comp_quarts = [sum(ptcomps{i}(first_quart)),...
                    sum(ptcomps{i}(second_quart)),...
                    sum(ptcomps{i}(third_quart)),...
                    sum(ptcomps{i}(fourth_quart))];
    comp_probs = comp_quarts./num_quarts;
 
    BetaInv84 = betainv(0.84,comp_quarts+1,num_quarts-comp_quarts+1);
    BetaInv16 = betainv(0.16,comp_quarts+1,num_quarts-comp_quarts+1);
    
    errorbar(meds,comp_probs,max(0,comp_probs-BetaInv16),max(0,BetaInv84-comp_probs),'k*','LineWidth',1);
    
    low_bins = [min(perfs{i}(first_quart)),...
                min(perfs{i}(second_quart)),...
                min(perfs{i}(third_quart)),...
                min(perfs{i}(fourth_quart))];
    high_bins = [max(perfs{i}(first_quart)),...
                max(perfs{i}(second_quart)),...
                max(perfs{i}(third_quart)),...
                max(perfs{i}(fourth_quart))];
            
    errorbar_x(meds,comp_probs,(meds-low_bins),(high_bins-meds),'k*');  
    text(0.01,0.45,['log-likelihood = $',num2str(probit_llhds(i),4),'$',10,...
        '$p$-value = $',num2str(probit_pvals(i),2),'$'],...
        'interpreter','latex','fontsize',24);
    set(gca,'FontSize',16);
    xlabel(['Perfusion Reduction (%)'],'FontSize',20);
    ylabel(['Probability of Complication (%)'],'FontSize',20);
     
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'rp_probit_resp_',grp_labels{i}],'-pdf');
        disp(['Saving ',fig_loc,'rp_probit_resp_',grp_labels{i},'.pdf...']);
    end
    
    %% hitogram of perf values for each dataset
    cur_fig=figure(i+1+2*length(groups)); clf reset; 
    set(cur_fig,'Position',ss_four2three);
   
%     
%     hist(perfs{i}(logical(ptcomps{i})),50);hold on;
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','r');
%     hist(perfs{i}(~logical(ptcomps{i})),50);
%     set(gca,'FontSize',18);
    
    [n1,x1] =hist(perfs{i}(~logical(ptcomps{i})),50);    
    [n2,x2]=hist(perfs{i}(logical(ptcomps{i})),50);

    
    bar(x1,n1,'hist');hold on;
    h=bar(x2,n2,'hist');hold off;
    set(h,'facecolor','r');
     set(gca,'FontSize',18);
     
    xlabel('Perfusion Reduction (%)','FontSize',20);
    ylabel('DVHs','FontSize',20);
    if do_print,
        set(cur_fig,'Color','w');
        export_fig(cur_fig,[fig_loc,'h_',grp_labels{i},'_perfs'],'-pdf');
        disp(['Saving ',fig_loc,'h_',grp_labels{i},'_perfs.pdf...']);
    end
    
    
end

