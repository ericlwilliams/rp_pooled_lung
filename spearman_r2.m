function spearman_r2
tic;
    %fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_EUD_meta.mat'; % RTOG
    %data
    fn='Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_RTOG_EUD_meta.mat'; % RTOG data'; % RTOG data
    load(fn,'CGmsk','CGcomb','CGnki','CGrtog');
     
    rtog_grp = CGrtog(1).mGrp;
    msk_grp = CGmsk(1).mGrp;
    nki_grp = CGnki(1).mGrp;
    comb_grp = CGcomb(1).mGrp;
    
    grps = {msk_grp nki_grp rtog_grp comb_grp};
    groups = {'MSK' 'NKI' 'RTOG' 'COMB'};
    
    x_axis = -1:0.1:1;
    
    %msk_euds = [grps{1}.mEUD]';
    %nki_euds = [grps{2}.mEUD]';
    %rtog_euds = [grps{3}.mEUD]';
    %comb_euds = [grps{4}.mEUD]';
    
    for i=1:length(grps),
        comps = [grps{i}.mFlgCensor]';
        comps = ~comps;
        euds = [grps{i}.mEUD]';
        rho = corr(euds,comps,'type','Spearman');
        Rs = rho.^2;
        
        figure(i);
        ylim([0 0.08]);
        % reverse order of Rs to display in log(a) instead of log(n)
        plot(x_axis,flipud(Rs),'*');
        xlabel('log10a');
        ylabel('Rs');
        grid;
    end
        %% Spearman's 
    %rtog_comps = [rtog_grp.mFlgCensor]';
    %rtog_comps = ~rtog_comps;
    %rtog_euds= [rtog_grp.mEUD]';
    %rtog_rho = corr(rtog_euds,rtog_comps,'type','Spearman');
    %rtog_Rs = rtog_rho.^2;
    %x_axis = -1:0.1:1;
    %ylim([0 0.08]);
    %plot(x_axis,rtog_Rs,'*');
    %xlabel('log10a');
    %ylabel('Rs');
    %grid;
    
    toc;
    
end