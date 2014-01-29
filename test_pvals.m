function test_pvals

fn = 'Z:\elw\MATLAB\original_data\MSK\rp_comb_best_euds.xlsx';

xlsread(fn);

inst = ans(:,1);
rp = ans(:,2);
eud = ans(:,3);
all = ones(length(eud),1);
%probit
[prb_b,~,prb_stats] = glmfit(eud,[rp all],'binomial','link','probit');
prb_pvals = prb_stats.p;
prb_p = prb_pvals(2);

yfit = glmval(prb_b,eud,'probit',prb_stats);
disp(['Probit p-value for log10a=-0.2 => ',num2str(prb_p)]);


[~,~,lgr_stats] = glmfit(eud,[rp all],'binomial','link','logit');
lgr_pvals = lgr_stats.p;
lgr_p = lgr_pvals(2);

disp(['Logit p-value for log10a=-0.2 => ',num2str(lgr_p)]);




end