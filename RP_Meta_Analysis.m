function RP_Meta_Analysis
% do meta analysis (fixed/random effects) for rp study
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
%             msk,   nki,   rtog,  umich
lkb_log10a = [-0.34, -1,  0.33, 0.09];
lkb_log10a_lcl = [-Inf, -Inf, -0.126, -0.199];
lkb_log10a_ucl = [0.137, 0.48, 0.827, 0.252];

lkb_log10n = -lkb_log10a;
lkb_log10n_lcl = -lkb_log10a_lcl;
lkb_log10n_ucl = -lkb_log10a_ucl;


disp(['a: ',num2str(10.^(lkb_log10a))]);
disp(['a lcl: ',num2str(10.^(lkb_log10a_lcl))]);
disp(['a ucl: ',num2str(10.^(lkb_log10a_ucl))]);
disp(' ');
disp(['n: ',num2str(10.^(lkb_log10n))]);
disp(['n lcl: ',num2str(10.^(lkb_log10n_ucl))]);
disp(['n ucl: ',num2str(10.^(lkb_log10n_lcl))]);

%return
% inst x param
%            a,    n,    td50, m
lkb_vals = [[0.457, 2.188, 14.9, 0.42];... %msk
            [0.1,  10,   17,   0.61];... % nki
            [2.138, 0.468, 45.9, 0.21];... % rtog
            [1.23, 0.81, 22.1, 0.16]]; % umich
        
        % errors with corrected log10(a/n) -> a/n
        % errors on n from invertign errors on a         
   lkb_lcl = [[0.1,  0.729, 7.9,  0.24];...
           [0.1,  0.331, 11.5, 0.41];...
           [0.748, 0.14894, 20.3, 0.12];...
           [0.6324, 0.55976, 16.1, 0.10]];
       
       lkb_ucl = [[1.37, 10,   42.0,  0.78];...
           [3.02, 10,   100.0, 1.09];...
           [6.7143, 1.33, 87.2,  0.41];...
           [1.787, 1.32, 29.1,  0.29]];
       
   
       
%        lkb_vals = [[0.46, 2.17, 14.9, 0.42];... %msk
%             [0.1,  10,   17,   0.61];... % nki
%             [2.14, 0.47, 45.9, 0.21];... % rtog
%             [1.23, 0.81, 22.1, 0.16]]; % umich
%         
% %         % errors on n from invertign errors on a         
%     lkb_ucl = [[1.37, 10,   42.0,  0.78];...
%            [3.02, 10,   100.0, 1.09];...
%            [6.72, 1.33, 87.2,  0.41];...
%            [1.79, 1.32, 29.1,  0.29]];
%        
%     lkb_lcl = [[0.1,  0.73, 7.9,  0.24];...
%            [0.1,  0.33, 11.5, 0.41];...
%            [0.75, 0.15, 20.3, 0.12];...
%            [0.76, 0.56, 16.1, 0.10]];
    

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
       
       
pld_vals = [0.69,1.45,19.2,0.38];
pld_lcl = [0.33,0.94,14.3,0.29];
pld_ucl = [1.06,3.03,26.2,0.51];

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



 




end

