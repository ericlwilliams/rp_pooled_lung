function RP_Anonymize

% fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat';
%  fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_med_EUD_fine_meta.mat';

eud_binning = 'fine'; % log10a binning
 
fns={['Z:/elw/MATLAB/meta_analy/meta_data/MSK_',eud_binning,'_EUD_meta.mat'],...
    ['Z:/elw/MATLAB/meta_analy/meta_data/NKI_',eud_binning,'_EUD_meta.mat'],...
    ['Z:/elw/MATLAB/meta_analy/meta_data/UMich_',eud_binning,'_EUD_meta.mat'],...
    ['Z:/elw/MATLAB/meta_analy/meta_data/RTOG_',eud_binning,'_EUD_meta.mat'],...
    ['Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta_comb.mat'],...
    ['Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_fine_EUD_fine_meta.mat']};

cohorts = {'MSK','NKI','UMICH','RTOG','MSK_NKI_UMich_RTOG_comb','MSK_NKI_UMich_RTOG'};

for i=1:length(fns)
    
    disp(['Anonymizing ',cohorts{i},'...']);
    clear CGobjs;
    load(fns{i});
    clear CGstrct;

    if isequal(cohorts{i},'NKI')
        CGobjs = CGobjs(1); % regions
    end
    if isequal(cohorts{i},'MSK_NKI_UMich_RTOG') ||...
            isequal(cohorts{i},'MSK_NKI_UMich_RTOG_comb')
        CGobjs = CGcomb;
        clear CGcomb;
    end
    
    CGgrp = CGobjs.mGrp;

    % encrypt mID 
   anon_mrns = cell(length(CGgrp),2);
   for j=1:length(CGgrp)
       
       
       tmp_id = CGgrp(j).mID;
        tmp_md5 = DataHash(tmp_id);
        anon_mrns{j,1} = tmp_id;
        anon_mrns{j,2} = tmp_md5;
       
        CGgrp(j).mID = tmp_md5;

    % Remove dates    
    if ~isempty([CGgrp(j).mDateComp]) ||...
       ~isempty([CGgrp(j).mDateStartTx]) ||...
       ~isempty([CGgrp(j).mDateEndTx]) ||...
       ~isempty([CGgrp(j).mDateLastFollowUp]) ||...
       ~isempty([CGgrp(j).mDateCensor]) ||...
       ~isempty([CGgrp(j).mDateComp]) ||...
       ~isempty([CGgrp(j).mDateRelapse]) ||...
       ~isempty([CGgrp(j).mDateBirth]) ||...
       ~isempty([CGgrp(j).mDateDeath]),
           
                CGgrp(j).mDateComp = [];
                CGgrp(j).mDateStartTx = [];
                CGgrp(j).mDateEndTx = [];
                CGgrp(j).mDateLastFollowUp = [];
                CGgrp(j).mDateCensor = [];
                CGgrp(j).mDateComp = [];
                CGgrp(j).mDateRelapse = [];
                CGgrp(j).mDateBirth = [];
                CGgrp(j).mDateDeath = [];
                
                
            end
   end

   [~,md5_inds] = sort({anon_mrns{:,2}});
   new_ids = mat2cell(md5_inds,1,ones(1,size(md5_inds,2)));
   anon_mrns(:,2) = new_ids';
   
   for k=1:length(CGgrp)
       CGgrp(k).mID = num2str(anon_mrns{k,2});
   end
   CGobjs.mGrp = CGgrp;
   xlswrite(['Z:/elw/MATLAB/meta_analy/meta_data/anon_data/',cohorts{i},'_mrn_encrypt.xlsx'],anon_mrns)
   anon_fn=['Z:/elw/MATLAB/meta_analy/meta_data/anon_data/',cohorts{i},'_',eud_binning,'_EUD_meta_anon.mat'];

   if isequal(cohorts{i},'MSK_NKI_UMich_RTOG') ||...
            isequal(cohorts{i},'MSK_NKI_UMich_RTOG_comb')
        CGcomb = CGobjs;
        clear CGobjs
        save(anon_fn,'CGcomb');
   else
        save(anon_fn,'CGobjs');
   end
end
end
    
