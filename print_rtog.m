function print_rtog

    load('Z:/elw/MATLAB/meta_analy/meta_data/RTOG_EUD_meta.mat');
    
    
    mNumPts = CGobjs(1).mNumInGrp;
    pt_grp = CGobjs(1).mGrp;
    
    tot_doses = zeros(mNumPts,1);
    for i=1:mNumPts,
        cur_pt = pt_grp(i);
        disp(['%%%%%%%%%%%%%%%']);
        disp(['Patient: ', cur_pt.mID])
        disp(['  Tx Date: ', datestr(cur_pt.mDateBaseline)]);
        
        %if cur_pt.mCompGrade,
            disp(['  Comp Grade: ', num2str(cur_pt.mCompGrade)]);
            disp(['  Comp Date: ', datestr(cur_pt.mDateComp)]);
        %end
        
    end
    
end