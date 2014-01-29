function qa_rtog
tic;
    %fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/rtog_hetero_EUD_meta.mat'; % RTOG
    %data
    fn='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_Hetero_EUD_meta_comp_lt6m.mat'; % RTOG data'; % RTOG data
    load(fn,'CGobjs');
     
    hetero_grp = CGobjs(1).mGrp;
    
    fn='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_Homo_EUD_meta_comp_lt6m.mat'; % RTOG data'; % RTOG data
    load(fn,'CGobjs');
    
    homo_grp = CGobjs(1).mGrp;
    
    x_axis = -1:0.1:1;
    
            %% Spearman's 

    
    hetero_comps = [hetero_grp.mFlgCensor]';
    hetero_comps = ~hetero_comps; 

    hetero_euds= [hetero_grp.mEUD]';
    hetero_rho = corr(hetero_euds,hetero_comps,'type','Spearman');

    homo_comps = [homo_grp.mFlgCensor]';
    homo_comps = ~homo_comps;  
    homo_euds= [homo_grp.mEUD]';
    homo_rho=corr(homo_euds,homo_comps,'type','Spearman');
        
    figure(1);
    hold on;
   
    plot(-x_axis,hetero_rho,'b*');
    plot(-x_axis,homo_rho,'r*');
   
    %ylim([0 0.3]); 
    xlabel('log10a');
    ylabel('Rs');
    grid;
    
  
    toc;
    
end