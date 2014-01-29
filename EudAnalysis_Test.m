function EudAnalysis_Test
tic; %close all;
fig_loc = 'Z:/elw/MATLAB/meta_analy/figures/latest/';
fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_meta.mat';
load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');

msk = CGmsk.mGrp;
msk_comps = ~[msk(:).mFlgCensor];
nki = CGnki.mGrp;
nki_comps = ~[nki(:).mFlgCensor];
umich = CGum.mGrp;
umich_comps = ~[umich(:).mFlgCensor];
rtog = CGrtog.mGrp;
rtog_comps = ~[rtog(:).mFlgCensor];

disp('n_pts & n_comps & Overall Incidence')
disp(['MSK: ',num2str(length(msk)),' & ',...
    num2str(length(msk(msk_comps))),' & ',...
    num2str(length(msk(msk_comps))/length(msk))]);
disp(['NKI: ',num2str(length(nki)),' & ',...
    num2str(length(nki(nki_comps))),' & ',...
    num2str(length(nki(nki_comps))/length(nki))]);
disp(['RTOG: ',num2str(length(rtog)),' & ',...
    num2str(length(rtog(rtog_comps))),' & ',...
    num2str(length(rtog(rtog_comps))/length(rtog))]);
disp(['UMich: ',num2str(length(umich)),' & ',...
    num2str(length(umich(umich_comps))),' & ',...
    num2str(length(umich(umich_comps))/length(umich))]);



disp([]);

end