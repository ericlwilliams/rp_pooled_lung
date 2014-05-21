function RP_Meta_Analysis
% do meta analysis (fixed/random effects) for rp study
% does forest plots
% saves meta analysis results to lkb_meta_parameters.mat


screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
do_one_sided_vars = false;
do_propigate_n_errors = false;

disp(['== Settings ==']);
disp(['do_print: ',num2str(do_print)]);
disp(['do_one_sided_vars: ',num2str(do_one_sided_vars)]);
disp(['do_propigate_n_errors: ',num2str(do_propigate_n_errors)]);
disp('==========');

fig_loc = 'Z:\elw\MATLAB\meta_analy\slides\figures\latest\';
lkb_loc = 'Z:/elw/MATLAB/meta_analy/meta_data/lkb_parameters.mat';
save_loc = 'Z:/elw/MATLAB/meta_analy/meta_data/lkb_meta_parameters.mat';

if isunix %on mac
    fig_loc = strrep(fig_loc,'Z:/elw/','/Users/elw/Documents/');
    lkb_loc = strrep(lkb_loc,'Z:/elw/','/Users/elw/Documents/');
    save_loc = strrep(save_loc,'Z:/elw/','/Users/elw/Documents/');
end

% lkb values from RP_Pooled_Results.m
load(lkb_loc,'lkb_vals','lkb_lcl','lkb_ucl');


pld_vals = lkb_vals(5,:);  
pld_lcl = lkb_lcl(5,:);  
pld_ucl = lkb_ucl(5,:);  

lkb_vals = lkb_vals(1:4,:);
lkb_lcl = lkb_lcl(1:4,:);
lkb_ucl = lkb_ucl(1:4,:);
  

if do_propigate_n_errors, % recalc errors on n from error propigation
    if ~do_one_sided_vars,
        disp(['TMP: only propigating n errors with one-sided variances']);
    end
    
    
    lkb_a = lkb_vals(:,1);
    
    
    lkb_a_lcl = lkb_lcl(:,1);
    lkb_a_ucl = lkb_ucl(:,1);

    % where lcl hits boundary (0.1), symmetrize from upper     
%     lkb_a_lcl(lkb_a_lcl==0.1) =...
%         lkb_a(lkb_a_lcl==0.1) - (lkb_a_ucl(lkb_a_lcl==0.1) - lkb_a(lkb_a_lcl==0.1));
    
    % ste
    lkb_a_lse = (lkb_a - lkb_a_lcl)./(1.96*lkb_a);
    lkb_a_use = (lkb_a_ucl - lkb_a)./(1.96*lkb_a);
    
    lkb_n = lkb_vals(:,2);
    lkb_n_lse = lkb_a_lse.*lkb_n;
    lkb_n_use = lkb_a_use.*lkb_n;
    
    lkb_n_ucl = lkb_n + (1.96*lkb_n_use);
    lkb_n_lcl = lkb_n - (1.96*lkb_n_lse);
    
end

if do_one_sided_vars,
    % choose which half of ci to use
    % done by eye (see elw_meta_maps_2014_04_28.pdf)
    % 0 - left CI, 1 - right CI
    lkb_sig_ci = [[1,0,1,0];...
                  [1,0,1,0];...
                  [0,1,0,1];...
                  [0,1,0,1]];
    lkb_sig = (lkb_vals-lkb_lcl).*(~lkb_sig_ci) + (lkb_ucl-lkb_vals).*lkb_sig_ci;          
    lkb_sig = lkb_sig./1.96; % to ste
else % two sided vars, use half full width of sigma for variances
    lkb_sig = abs(lkb_ucl-lkb_lcl)./(2*1.96);
end

lkb_vars = lkb_sig.^2;
lkb_weights = 1./lkb_vars;

% fixed effect
fe_vals = sum(lkb_vals.*lkb_weights)./sum(lkb_weights);
fe_vars = 1./sum(lkb_weights);

fe_ucl = fe_vals + 1.96*sqrt(fe_vars);
fe_lcl = fe_vals - 1.96*sqrt(fe_vars);


% random effects
re_df = 4-1; %(number of studies) -1

re_q = sum(lkb_weights.*(lkb_vals.^2)) - ...
        (sum(lkb_weights.*lkb_vals).^2)./(sum(lkb_weights));

re_q_pval = 1-chi2cdf(re_q,re_df);
re_i2 = (re_q - re_df)./re_q;
re_i2(re_i2 < 0) = 0;

re_c = sum(lkb_weights) - (sum(lkb_weights.^2)./sum(lkb_weights));
    


re_tau_sq = (re_q - re_df)./re_c;
re_tau_sq(re_q < re_df) = 0;


re_vars = lkb_vars+[re_tau_sq;re_tau_sq;re_tau_sq;re_tau_sq];
re_weights = 1./re_vars;

re_vals = sum(lkb_vals.*re_weights)./sum(re_weights);
re_vars = 1./sum(re_weights);

re_ucl = re_vals + 1.96*sqrt(re_vars);
re_lcl = re_vals - 1.96*sqrt(re_vars);

%%%%%%%%%%%%%%%%%%%%%%%%
%% forest plots
%%%%%%%%%%%%%%%%%%%%%%%%
texts = cell(7,1);
texts{1} = 'MSK';
texts{2} = 'NKI';
texts{3} = 'RTOG';
texts{4} = 'UMICH';
texts{5} = 'FEM';
texts{6} = 'REM';
texts{7} = 'Pooled';

lkb_params = {'a','n','TD50','m'};

disp(['Meta analysis results (used in RP_Meta_Maps.m']);
disp(['Parameter |',' Fixed-Effects | ',' Random-Effects']);
 
mkr_size = repmat(5,length(texts),1);
for j=1:length(pld_vals)
    lkb_meta_vals = [lkb_vals(:,j);fe_vals(j);re_vals(j);pld_vals(j)];
    lkb_meta_lcl = [lkb_lcl(:,j);fe_lcl(j);re_lcl(j);pld_lcl(j)];
    lkb_meta_ucl = [lkb_ucl(:,j);fe_ucl(j);re_ucl(j);pld_ucl(j)];
   
    cur_fig=figure(j); clf reset; hold on;
    set(cur_fig,'Position',ss_four2three);
    
    forest(lkb_params{j},flipud(texts),flipud(lkb_meta_vals),...
        flipud(lkb_meta_lcl),flipud(lkb_meta_ucl),flipud(mkr_size),...
        re_q(j),re_q_pval(j),re_tau_sq(j),re_i2(j));
        
    
    if do_print,
        set(cur_fig,'Color','w');
%         set(cur_fig,'Renderer','opengl');
        export_fig(cur_fig,[fig_loc,'forest_',lkb_params{j}],'-pdf');
%         disp(['Saving ',fig_loc,'forest_',lkb_params{j},'.pdf...']);
    end
%% Printing results
disp(' ');
disp([lkb_params{j},', ',...
    'Q = ',num2str(re_q(j),4),', ',...
     'p = ',num2str(re_q_pval(j),4),', ',...
     'T2 = ',num2str(re_tau_sq(j),4),', ',...
     'I2 = ',num2str(re_i2(j),4)]);
    disp(['fe_',lkb_params{j},' = [',...
    num2str(fe_vals(j),4),',',num2str(fe_lcl(j),4),',',num2str(fe_ucl(j),4),']']);
disp(['re_',lkb_params{j},' = [',...
    num2str(re_vals(j),4),',',num2str(re_lcl(j),4),',',num2str(re_ucl(j),4),']']);

end



 disp(['Saving meta-analysis lkb parameters to: ',save_loc,'...']);
 save(save_loc,'fe_vals','fe_lcl','fe_ucl',...
               're_vals','re_lcl','re_ucl');




end

