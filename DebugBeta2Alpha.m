function DebugBeta2Alpha
tic;



load('MSK_NKI_RTOG_ntd_a2b_Hetero_EUD_meta_lt6m.mat');
rtog_ntd = CGrtog;
load('MSK_NKI_RTOG_bed_a2b_Hetero_EUD_meta_lt6m.mat');
rtog_bed = CGrtog;

disp(['']);

load('RTOG_Hetero_EUD_meta_comp_lt6m_phys_a2b.mat');
phys = CGobjs.mGrp;
phys_euds = [phys.mEUD];
phys_dosebins_lq = [phys.mDoseBins_LQ];
phys_dosebins_org = [phys.mDoseBins_org];


load('RTOG_Hetero_EUD_meta_comp_lt6m_ntd_a2b.mat');
ntd = CGobjs.mGrp;
ntd_euds = [ntd.mEUD];
ntd_dosebins_lq = [ntd.mDoseBins_LQ];
ntd_dosebins_org = [ntd.mDoseBins_org];
ntd_mld = ntd_euds(11,:);

load('RTOG_Hetero_EUD_meta_comp_lt6m_bed_a2b.mat');
bed = CGobjs.mGrp;
bed_euds = [bed.mEUD];
bed_dosebins_lq = [bed.mDoseBins_LQ];
bed_dosebins_org = [bed.mDoseBins_org];
bed_mld = bed_euds(11,:);

ntd_x=1:.05:21;
figure(1);
hist(ntd_mld,ntd_x);
figure(2);
bed_x=1:.05:35;
hist(bed_mld,bed_x);

disp(['']);


toc;
end
