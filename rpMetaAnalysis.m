function rpMetaAnalysis
% do meta analysis (fixed/random effects) for rp study
screen_size=get(0,'ScreenSize');
ss_four2three = [0 0 screen_size(3)/2 (screen_size(4)/2)*(4/3)];

do_print = true;
fig_loc = 'Z:\elw\MATLAB\meta_analy\slides\figures\latest\';
% inst x param
%            a,    n,    td50, m
lkb_vals = [[0.46, 2.17, 14.9, 0.42];... %msk
            [0.1,  10,   17,   0.61];... % nki
            [2.14, 0.47, 45.9, 0.21];... % rtog
            [1.23, 0.81, 22.1, 0.16]]; % umich
     
lkb_ucl = [[1.37, 10,   42.0,  0.78];...
           [3.02, 10,   100.0, 1.09];...
           [6.72, 1.33, 87.2,  0.41];...
           [1.79, 1.32, 29.1,  0.29]];
       
lkb_lcl = [[0.1,  0.73, 7.9,  0.24];...
           [0.1,  0.33, 11.5, 0.41];...
           [0.75, 0.15, 20.3, 0.12];...
           [0.76, 0.56, 16.1, 0.10]];

pld_vals = [0.69,1.45,19.2,0.38];
pld_lcl = [0.33,0.94,14.3,0.29];
pld_ucl = [1.06,3.03,26.2,0.51];
       
% choose which half of ci to use
% done by eye (see elw_meta_maps_2014_04_28.pdf)
% 0 - left CI, 1 - right CI
lkb_sig_ci = [[1,0,1,0];...
              [1,0,1,0];...
              [0,1,0,1];...
              [0,1,0,1]];
lkb_sig = (lkb_vals-lkb_lcl).*(~lkb_sig_ci) + (lkb_ucl-lkb_vals).*lkb_sig_ci;          
lkb_sig = lkb_sig./1.96; % to ste

% use half full width of sigma for variance
% lkb_sig = abs(lkb_ucl-lkb_lcl)./(2*1.96);
lkb_vars = lkb_sig.^2;

lkb_weights = 1./lkb_vars;

% fixed effect
fe_vals = sum(lkb_vals.*lkb_weights)./sum(lkb_weights);
fe_vars = 1./sum(lkb_weights);

fe_ucl = fe_vals + 1.96*sqrt(fe_vars);
fe_lcl = fe_vals - 1.96*sqrt(fe_vars);


% random effects

re_q = sum(lkb_weights.*(lkb_vals.^2)) - ...
        (sum(lkb_weights.*lkb_vals).^2)./(sum(lkb_weights));

re_c = sum(lkb_weights) - (sum(lkb_weights.^2)./sum(lkb_weights));
    
re_df = 3; %(number of studies) -1

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

mkr_size = repmat(5,length(texts),1);
% LKB a
for j=1:length(pld_vals)
    lkb_meta_vals = [lkb_vals(:,j);fe_vals(j);re_vals(j);pld_vals(j)];
    lkb_meta_lcl = [lkb_lcl(:,j);fe_lcl(j);re_lcl(j);pld_lcl(j)];
    lkb_meta_ucl = [lkb_ucl(:,j);fe_ucl(j);re_ucl(j);pld_ucl(j)];
   
    cur_fig=figure(j); clf reset; hold on;
    set(cur_fig,'Position',ss_four2three);
    
    forest(lkb_params{j},flipud(texts),flipud(lkb_meta_vals),...
        flipud(lkb_meta_lcl),flipud(lkb_meta_ucl),flipud(mkr_size));
    
    if do_print,
        set(cur_fig,'Color','w');
%         set(cur_fig,'Renderer','opengl');
        export_fig(cur_fig,[fig_loc,'forest_',lkb_params{j}],'-pdf');
        disp(['Saving ',fig_loc,'forest_',lkb_params{j},'.pdf...']);
    end
    
    end



end

