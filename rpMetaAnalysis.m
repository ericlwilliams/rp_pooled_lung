function rpMetaAnalysis
% do meta analysis (fixed/random effects) for rp study

% msk,nki,rtog,umich,comb
lkb_a = [0.46,0.1,2.14,1.23, 0.69];
lkb_a_ll = [0.1,0.1,0.75,0.76,0.33];
lkb_a_ul = [1.37,3.02,6.72,1.79,1.06];

lkb_n = [2.17,10,0.47,0.81,1.45];
lkb_n_ll = [0.73,0.33,0.15,0.56,0.94];
lkb_n_ul = [10,10,1.33,1.32,3.03];

lkb_td50 = [14.9,17,45.9,22.1,19.2];
lkb_td50_ll = [7.9,11.5,20.3,16.1,14.3];
lkb_td50_ul = [42.0,100,87.2,29.2,26.2];

lkb_m = [0.42,0.61,0.21,0.16,0.38];
lkb_m_ll = [0.24,0.41,0.12,0.10,0.29];
lkb_m_ul = [0.78,1.09,0.41,0.29,0.51];

% 
% pooled_vals = [1.45,0.69,19.2,0.38);
% pooled_ll = [0.94,0.33,14.3,0.29);
% pooled_ul = [3.03,1.06,26.2,0.51);

% sms_names = {'MSK','NKI','RTOG','UMICH','COMB'};
sms_names = {'MSK','NKI','RTOG','UMICH'};
sms = repmat(StatModel(),length(sms_names));

for i=1:length(sms_names)
    
    sms(i) = sms(i).SetModelName(sms_names{i});
    sms(i) = sms(i).SetModelParamIDs({'a','n','td50','m'});
    sms(i) = sms(i).SetModelParamVals([lkb_a(i),lkb_n(i),lkb_td50(i),lkb_m(i)]);
    sms(i) = sms(i).SetModelUncertainties(...
        [lkb_a_ll(i),lkb_n_ll(i),lkb_td50_ll(i),lkb_m_ll(i)],...
        [lkb_a_ul(i),lkb_n_ul(i),lkb_td50_ul(i),lkb_m_ul(i)]);
end



% 
% sm = sm.SetModelName('LKB');
% sm = sm.SetModelParamIDs({'a','n','td50','m'});
% sm = sm.SetModelParamVals([1, 2, 3, 4]);
% % 95% CL
% sm = sm.SetModelUncertainties(...
%     [0.5,1.5,2.5,3.5],...
%     [1.5,2.5,3.5,4.5]);
% 
% 


end