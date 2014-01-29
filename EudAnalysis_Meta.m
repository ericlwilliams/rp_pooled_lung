function EudAnalysis_Meta
tic; %close all;

comb_only = false;
ppm_only = true;

binning = 'fine'; % TD50 and m binning
eud_binning = 'med';
% load results
fn = ['Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_',eud_binning,'_EUD_',binning,'_meta.mat'];
%fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_RTOG_Homo_EUD_meta.mat';

% reload MSK/NKI data (outputs from Pt_MSK/NKI, combines them
% into
% EUD_regional_MSK_NKI.mat, and saves analysis data as calculated below
if 1
    
   
    if ppm_only
        fn_msk=['Z:/elw/MATLAB/meta_analy/meta_data/MSK_',eud_binning,'_EUD_meta_ppm.mat']; % MSK data - ONLY coarse log(n) available
        fn_nki=['Z:/elw/MATLAB/meta_analy/meta_data/NKI_',eud_binning,'_EUD_meta_ppm.mat']; % NKI data
        fn_um=['Z:/elw/MATLAB/meta_analy/meta_data/UMich_',eud_binning,'_EUD_meta_ppm.mat']; % NKI data
        fn_rtog=['Z:/elw/MATLAB/meta_analy/meta_data/RTOG_',eud_binning,'_EUD_meta_ppm.mat'];
    else
        fn_msk=['Z:/elw/MATLAB/meta_analy/meta_data/MSK_',eud_binning,'_EUD_meta.mat']; % MSK data - ONLY coarse log(n) available
        fn_nki=['Z:/elw/MATLAB/meta_analy/meta_data/NKI_',eud_binning,'_EUD_meta.mat']; % NKI data
        fn_um=['Z:/elw/MATLAB/meta_analy/meta_data/UMich_',eud_binning,'_EUD_meta.mat']; % NKI data
        fn_rtog=['Z:/elw/MATLAB/meta_analy/meta_data/RTOG_',eud_binning,'_EUD_meta.mat'];
    end
    %data
    %fn_rtog='Z:/elw/MATLAB/meta_analy/meta_data/RTOG_Homo_EUD_meta.mat'; % RTOG data
    
    load(fn_rtog,'CGobjs'); CGrtog = CGobjs(1);
    load(fn_um,'CGobjs'); CGum = CGobjs;
    load(fn_nki,'CGobjs'); CGnki = CGobjs(1);
    load(fn_msk,'CGobjs'); CGmsk = CGobjs(1);
    
    
    CGcomb = CGrtog;
    CGcomb = CGcomb.fAddPatient(CGum.mGrp);
    CGcomb = CGcomb.fAddPatient(CGnki.mGrp);
    CGcomb = CGcomb.fAddPatient(CGmsk.mGrp);
    if ppm_only
        fn = ['Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_',eud_binning,'_EUD_',binning,'_meta_ppm.mat'];
    end
    save(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');
end

if comb_only
load(fn,'CGcomb');
else
load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');
end

% generate EUD atlases
if 0 && ~ppm_only
    eudstep = 0.5; % step size in atlas

    disp(['MSK ATLAS calculation']);
    CGmsk.mStepDose = eudstep;
    CGmsk = CGmsk.fCalculateEUDBins();
    CGmsk = CGmsk.fCrudeAtlas_EUD();
    
    disp(['NKI ATLAS calculation']);
    CGnki.mStepDose = eudstep;
    CGnki = CGnki.fCalculateEUDBins();
    CGnki = CGnki.fCrudeAtlas_EUD();

    disp(['UMich ATLAS calculation']);
    CGum.mStepDose = eudstep;
    CGum = CGum.fCalculateEUDBins();
    CGum = CGum.fCrudeAtlas_EUD();
    
    disp(['RTOG ATLAS calculation']);
    CGrtog.mStepDose = eudstep;
    CGrtog = CGrtog.fCalculateEUDBins();
    CGrtog = CGrtog.fCrudeAtlas_EUD();
    
    disp(['COMB ATLAS calculation']);
    CGcomb.mStepDose = eudstep;
    CGcomb = CGcomb.fCalculateEUDBins();
    CGcomb = CGcomb.fCrudeAtlas_EUD();
    
end

% Atlas analysis - probability of having >=20% RP rate & low 68% confidence
if 0 && ~ppm_only
    % probability
    disp(['MSK beta cumulative']);
    CGmsk = CGmsk.fBetaCumulativeProbability_EUD();
    disp(['NKI beta cumulative']);
    CGnki = CGnki.fBetaCumulativeProbability_EUD();
    disp(['UMich beta cumulative']);
    CGum = CGum.fBetaCumulativeProbability_EUD();
    disp(['RTOG beta cumulative']);
    CGrtog = CGrtog.fBetaCumulativeProbability_EUD();
    disp(['COMB beta cumulative']);
    CGcomb = CGcomb.fBetaCumulativeProbability_EUD();
    
    % confidence
    disp(['MSK beta inverse']);
    CGmsk = CGmsk.fBetaInverseProbability_EUD();
    disp(['NKI beta inverse']);
    CGnki = CGnki.fBetaInverseProbability_EUD();
    disp(['UMich beta inverse']);
    CGum = CGum.fBetaInverseProbability_EUD();
    disp(['RTOG beta inverse']);
    CGrtog = CGrtog.fBetaInverseProbability_EUD();
    disp(['COMB beta inverse']);
    CGcomb = CGcomb.fBetaInverseProbability_EUD();
    
end

% logistic regression analysis using matlab functions
if 1 && ~ppm_only
    
    if comb_only
    disp(['COMB Logistic Analysis']);    
    CGcomb = CGcomb.fLogisticRegressionExact_EUD();
        
        
    else
        
    disp(['MSK Logistic Analysis']);
    CGmsk = CGmsk.fLogisticRegressionExact_EUD();
    disp(['NKI Logistic Analysis']);
    CGnki = CGnki.fLogisticRegressionExact_EUD();
    disp(['UMich Logistic Analysis']);
    CGum = CGum.fLogisticRegressionExact_EUD();
    disp(['RTOG Logistic Analysis']);    
    CGrtog = CGrtog.fLogisticRegressionExact_EUD();
    disp(['COMB Logistic Analysis']);    
    CGcomb = CGcomb.fLogisticRegressionExact_EUD();
    end
end
% logistic regression grid analysis
if 0 && ~ppm_only
    LogisticBetaRange = {(-10:0.1:10)'; (-1:0.01:1)'};
   
    CGmsk.mLogisticRegressionGridBetaRange = LogisticBetaRange;
    CGnki.mLogisticRegressionGridBetaRange = LogisticBetaRange;
    CGum.mLogisticRegressionGridBetaRange = LogisticBetaRange;
    CGrtog.mLogisticRegressionGridBetaRange = LogisticBetaRange;
    CGcomb.mLogisticRegressionGridBetaRange = LogisticBetaRange;
   
    disp(['MSK Logistic Grid']);            
    CGmsk = CGmsk.fLogisticRegressionGridExact_EUD();
    disp(['NKI Logistic Grid']);    
    CGnki = CGnki.fLogisticRegressionGridExact_EUD();
    disp(['UMich Logistic Grid']);    
    CGum = CGum.fLogisticRegressionGridExact_EUD();
    disp(['RTOG Logistic Grid']);    
    CGrtog = CGrtog.fLogisticRegressionGridExact_EUD();
    disp(['COMB Logistic Grid']);    
    CGcomb = CGcomb.fLogisticRegressionGridExact_EUD();
end

% logistic regression grid analysis - goodness of fit
if 0 && ~ppm_only
    disp(['MSK Logistic GoF']);    
    CGmsk = CGmsk.fLogisticRegressionGoodnessOfFitSimulationExact_EUD();
    disp(['NKI Logistic GoF']);    
    CGnki = CGnki.fLogisticRegressionGoodnessOfFitSimulationExact_EUD();
    disp(['UMich Logistic GoF']);    
    CGum = CGum.fLogisticRegressionGoodnessOfFitSimulationExact_EUD();
    disp(['RTOG Logistic GoF']);    
    CGrtog = CGrtog.fLogisticRegressionGoodnessOfFitSimulationExact_EUD();
    disp(['COMB Logistic GoF']);    
    CGcomb = CGcomb.fLogisticRegressionGoodnessOfFitSimulationExact_EUD();

end

% Lyman  analysis
if 0 && ~ppm_only
 
    disp(['RTOG Lyman Analysis']);
    CGrtog = CGrtog.fLymanAnalysisGridExact_EUD(binning);
    
    disp(['MSK Lyman Analysis']);
    CGmsk = CGmsk.fLymanAnalysisGridExact_EUD(binning);

    disp(['NKI Lyman Analysis']);
    CGnki = CGnki.fLymanAnalysisGridExact_EUD(binning);
    
    disp(['UMich Lyman Analysis']);
    CGum = CGum.fLymanAnalysisGridExact_EUD(binning);
    
    disp(['COMB Lyman Analysis']);
    CGcomb = CGcomb.fLymanAnalysisGridExact_EUD(binning);
    
end

% Lyman  analysis - goodness of fit
if 0 && ~ppm_only
        
    disp(['MSK Lyman GoF']);
    CGmsk = CGmsk.fLymanGoodnessOfFitSimulationExact_EUD();
    disp(['NKI Lyman GoF']);
    CGnki = CGnki.fLymanGoodnessOfFitSimulationExact_EUD();
    disp(['UMich Lyman GoF']);
    CGum = CGum.fLymanGoodnessOfFitSimulationExact_EUD();
    disp(['RTOG Lyman GoF']);
    CGrtog = CGrtog.fLymanGoodnessOfFitSimulationExact_EUD();
    disp(['COMB Lyman GoF']);
    CGcomb = CGcomb.fLymanGoodnessOfFitSimulationExact_EUD();
    
end

% save
if comb_only
    fn = ['Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_',eud_binning,'_EUD_',binning,'_meta_comb.mat'];
    save(fn,'CGcomb');
else
    if ppm_only
        fn = ['Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_',eud_binning,'_EUD_',binning,'_meta_ppm.mat'];
    end
    save(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');
end

toc;