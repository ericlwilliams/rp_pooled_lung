function rtog_ott

   fn='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_EUD_meta.mat'; % RTOG data'; % RTOG data
   load(fn);
    
   CGrtog = CGobjs;
    rtog_grp = CGrtog.mGrp;
    rtog_start_dates = [rtog_grp.mDateStartTx];
    rtog_end_dates = [rtog_grp.mDateEndTx];
    rtog_ott = rtog_end_dates-rtog_start_dates;
    
    boxplot(rtog_ott);
    xlabel('RTOG','FontSize',14);
    ylabel('OTT (days)','FontSize',14);
    text(1.1,50,['Median OTT: ',num2str(median(rtog_ott))],...
        'FontSize',14);
    ylim([0 100]);
end
