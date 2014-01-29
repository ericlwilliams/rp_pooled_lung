function load_rtog

to_save=false;

%% load_rtog
% this function opens RTOG-9311 whole lung data,
% applies alpha/beta = 3 Normalized Tissue Dose correction with 2 Gy
% fractions, and calculates GED values

tic;

    
eud_binning='fine';

%% Load RTOG data
%   classDataFromXls with properties
%       mXlsRaw: raw data of xls sheet
%        mLabel: column name
%        mData: data in the column/row
%        mFlg: flags of which column/row has the data
%
% CGobjs: (Complication Group objects) class contains all analysis info and
% data (atlases, proability results, lyman model parameters, etc.) for all
% patients (classOutcomeIndividual)
%   class: classOutComeAnalysis() with properties
%        mGrp: patient group of classOutcomeIndividual objects
%        mNumInGrp: number of patients in the group
%        mLymanN: range of n values
%

%missing_lung_dvhs = {'JACKSON_9311_092' 'JACKSON_9311_163' 'JACKSON_9311_166'};
missing_lung_dvhs = {'JACKSON_9311_092','JACKSON_9311_163',...
    'JACKSON_9311_006','JACKSON_9311_136' };

% currently excluding DVHs with < 1% change between Lung-PTV -> Lung-GTV
%excluded_dvhs = {'JACKSON_9311_031', 'JACKSON_9311_105',...
%    'JACKSON_9311_104', 'JACKSON_9311_121'};
 excluded_dvhs = {'JACKSON_9311_014'
    'JACKSON_9311_023'
    'JACKSON_9311_029'
    'JACKSON_9311_030'
    'JACKSON_9311_031'
    'JACKSON_9311_042'
    'JACKSON_9311_044'
    'JACKSON_9311_046'
    'JACKSON_9311_056'
    'JACKSON_9311_062'
    'JACKSON_9311_066'
    'JACKSON_9311_068'
    'JACKSON_9311_074'
    'JACKSON_9311_076'
    'JACKSON_9311_079'
    'JACKSON_9311_080'
    'JACKSON_9311_083'
    'JACKSON_9311_088'
    'JACKSON_9311_091'
    'JACKSON_9311_099'
    'JACKSON_9311_104'
    'JACKSON_9311_105'
    'JACKSON_9311_107'
    'JACKSON_9311_108'
    'JACKSON_9311_115'
    'JACKSON_9311_119'
    'JACKSON_9311_120'
    'JACKSON_9311_121'
    'JACKSON_9311_122'
    'JACKSON_9311_129'
    'JACKSON_9311_141'
    'JACKSON_9311_142'
    'JACKSON_9311_147'
    'JACKSON_9311_152'
    'JACKSON_9311_155'
    'JACKSON_9311_156'
    'JACKSON_9311_160'
    'JACKSON_9311_172'};

alpha2beta = 3;
a2b_corr = 'NTD';

PtInfo = classDataFromXls();
Regions = {'Whole'}';
CGobjs = classOutcomeAnalysis();


if isequal(eud_binning,'fine')
    CGobjs.mLymanN = 10.^(-1:0.01:1)';
elseif isequal(eud_binning,'med')
    CGobjs.mLymanN = 10.^(-1:0.02:1)';
elseif isequal(eud_binning,'crs')
    CGobjs.mLymanN = 10.^(-1:0.1:1)';
else %default coarse
    CGobjs.mLymanN = 10.^(-1:0.1:1)';
end


if ~strcmp(a2b_corr,'PHYS')
    CGobjs.mBeta2Alpha = [1./alpha2beta];
end

mNumRegions =1;% one analysis set for each lung region
CGobjs = repmat(CGobjs,[mNumRegions,1]);

%% read xls file
fn = 'Z:/elw/MATLAB/original_data/RTOG9311/2012_12_21/9311 Jackson data 20DEC2012.xlsx';
%fn = 'Z:/elw/MATLAB/original_data/RTOG9311/test.xlsx';
xlsPatientInfo=xlsFileRead(fn);
%fn ='Z:/elw/MATLAB/original_data/RTOG9311/2012_12_21/LUNG_TOTAL_Hetero_dvh - CORRECTED - 2012-12-20.xlsx';
fn ='Z:/elw/MATLAB/original_data/RTOG9311/2013_01_25/RTOG9311_LUNG-GTV_HETERO_DVH.xlsx';

xlsWholeLungInfo=xlsFileRead(fn);

PtInfo.mXlsRaw = xlsPatientInfo(1).xlsRaw; % load data (only 1 sheet)

%% Get number of patients and patient ids
PtInfo.mLabel = 'Patient ID'; % find the pt Ids
%# for RTOG data, patients are repeated for complications
PtInfo = PtInfo.fExtractColData();

flgPtId = PtInfo.mFlg; % pt with ID number are flaged
mFullPtIds = PtInfo.mData;
mFullPtIds = mFullPtIds(flgPtId);

%# Remove duplicate mPtIDs (from multiple lung toxicities
mPtIds = unique(mFullPtIds);
%# Number of patients
mNumPts = length(mPtIds);

%% Get censor info, start, end, and complication (>=3) date
PtInfo.mLabel = 'Lung Toxicity Grade'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullToxGrades = PtInfo.mData;
mFullToxGrades = mFullToxGrades(2:end);%Remove 'Lung Toxicity Grade'

PtInfo.mLabel = 'Date of Lung Toxicity'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullToxDates = PtInfo.mData;
mFullToxDates = mFullToxDates(2:end);%Remove 'Date of Lung Toxicity'

PtInfo.mLabel = 'Date RT Started'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullStartDates = PtInfo.mData;
mFullStartDates = mFullStartDates(2:end);%Remove 'Date of Lung Toxicity'

PtInfo.mLabel = 'Date RT Ended'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullEndDates = PtInfo.mData;
mFullEndDates = mFullEndDates(2:end);%Remove 'Date of Lung Toxicity'

PtInfo.mLabel = 'Date of Last Evaluation'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullLastFollowUp = PtInfo.mData;
mFullLastFollowUp = mFullLastFollowUp(2:end);

PtInfo.mLabel = 'Fractions'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullNumFractions = PtInfo.mData;
mFullNumFractions = mFullNumFractions(2:end);

PtInfo.mLabel = 'Prescribed Dose'; % find the pt Ids
PtInfo = PtInfo.fExtractColData();
mFullTxLevels = PtInfo.mData;
mFullTxLevels = mFullTxLevels(2:end);

mFullTxLevels(strcmp(mFullTxLevels,'Dose level1')) = {70.9};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level2')) = {77.4};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level3')) = {83.8};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level4')) = {90.3};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level5')) = {70.9};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level6')) = {77.4};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level7')) = {83.8};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level8')) = {64.5};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level9')) = {70.9};
mFullTxLevels(strcmp(mFullTxLevels,'Dose level10')) = {77.4};






%% Load DVH info
%# Change date in PtInfo to load DVHs
PtInfo.mXlsRaw = xlsWholeLungInfo(1).xlsRaw; % load data (only 1 sheet)
PtInfo.mLabel = 'DOSE';
PtInfo = PtInfo.fExtractRowData();
flgdose = PtInfo.mFlg;
dosebins = [cell2mat(PtInfo.mData(flgdose))];
dosebins = [dosebins;dosebins(end)+0.1];%add one more bin, corresponding to vol=0

mNumDoseBins = length(dosebins);


CIobjs = classOutcomeIndividual();


%CIobjs = repmat(CIobjs,mNumPts-length(missing_lung_dvhs),1); %
%mNumPts-1 - hack to remove pt92 and pt163
CIobjs = repmat(CIobjs,mNumPts,1); % mNumPts-1 - hack to remove pt92 and pt163


%# get maximum toxicity for each patient
omitted_etc=0;
omitted_dvh=0;
excluded_dvh=0;
omitted_lt6m=0;


n_comp=0;

cur_i=1;
%tmp
for i=1:mNumPts,
%      %tmp
%     if i==55
%         disp(mPtIds(i));
%     end
%for i=1:2,
    %# If DVH is missing (set by missing_lung_dvhs) skip
    dvh_missing = false;
    for j=1:length(missing_lung_dvhs),
        if strcmp(mPtIds(i),missing_lung_dvhs(j)),
            dvh_missing=true;
            
            break;
        end
    end
    if dvh_missing,
        omitted_dvh=omitted_dvh+1;
        continue;
    end
    
    dvh_excluded = false;
    for j=1:length(excluded_dvhs),
        if strcmp(mPtIds(i),excluded_dvhs(j)),
            dvh_excluded=true;
            break;
        end
    end
    if dvh_excluded,
        excluded_dvh=excluded_dvh+1;
        continue;
    end
    
    
    cur_pt = mPtIds(i);
    
   
    cur_inds = strmatch(cur_pt,mFullPtIds);
    
    %find max toxicity w/date
    all_tox = [mFullToxGrades{cur_inds}];
    all_tox_dates = [mFullToxDates(cur_inds)];
    all_end_dates = [mFullEndDates(cur_inds)];
    
    [max_tox, tox_ind] = max(all_tox);
    
    comp_grade=0;
    comp_ind=-1;
    if ~isnan(all_tox),
        for k=1:length(all_tox),
            days_until_tox =...
                (datenum(all_tox_dates(k))-datenum(all_end_dates(k)));
            if days_until_tox>=6*30,% more than 6 months for given toxicity
                break;
            else
                if all_tox(k)>comp_grade,
                    comp_grade=all_tox(k);
                    comp_ind=k;
                end
            end
        end
    else
        max_tox=0;
    end
    
    cur_max_tox_date = mFullToxDates(cur_inds(tox_ind));
    
    cur_start_date = mFullStartDates(cur_inds(tox_ind));
    cur_end_date = mFullEndDates(cur_inds(tox_ind));
    cur_followup_date = mFullLastFollowUp(cur_inds(tox_ind));
    cur_num_fractions = mFullNumFractions{cur_inds(tox_ind)};
    cur_tx = mFullTxLevels{cur_inds(tox_ind)};
    % pt 121 no start/end date
    if cur_start_date{1}=='.', % pt 121
        omitted_etc=omitted_etc+1;
        continue;
        %cur_start_date{1}='0/0/0';
    end
    
    if cur_end_date{1}=='.', % pt 121
        omitted_etc=omitted_etc+1;
        continue;
        %cur_end_date{1}='0/0/0';
    end
    
    lt6m_followup=false;
%if (datenum(cur_followup_date(1))-datenum(cur_start_date(1)))<6*30,
    if (datenum(cur_followup_date(1))-datenum(cur_end_date(1)))<6*30,
         lt6m_followup=true;
    end
    
    % Less than 6 months follow up, no complication >= 3Grade
    if lt6m_followup && comp_grade<3,
        omitted_lt6m=omitted_lt6m+1;
        continue;
    end
    
    
    if max_tox==0, % no complication recorded at all
        
        CIobjs(cur_i).mCompGrade = comp_grade;
        CIobjs(cur_i).mMaxToxGrade= max_tox;
        CIobjs(cur_i).mFlgCensor = true; %no complicaton - true
        CIobjs(cur_i).mDateCensor = datenum(cur_followup_date);

        CIobjs(cur_i).mDateComp = 0;
        
    else
        
        if comp_ind<0, % no complication recorded w/in 6 months
            
            CIobjs(cur_i).mCompGrade = comp_grade;%0
            CIobjs(cur_i).mMaxToxGrade= max_tox;
            CIobjs(cur_i).mFlgCensor = true; %no complicaton - true
            CIobjs(cur_i).mDateCensor = datenum(cur_followup_date); % date of maximum toxicity
            CIobjs(cur_i).mDateComp = datenum(cur_max_tox_date); % no complication w/in 6 months
            
        elseif  comp_grade < 3, %< 3 Grade within 6 months
            
            CIobjs(cur_i).mCompGrade = comp_grade;
            CIobjs(cur_i).mMaxToxGrade= max_tox;
            CIobjs(cur_i).mFlgCensor = true; % complication w/in 6 m
            CIobjs(cur_i).mDateCensor = datenum(mFullToxDates(cur_inds(comp_ind))); % date of maximum toxicity
            CIobjs(cur_i).mDateComp = datenum(cur_max_tox_date); % complication w/in 6 m
            
            
        else % >= 3Grade within 6 months
            n_comp=n_comp+1;
            
            CIobjs(cur_i).mCompGrade = comp_grade;
            CIobjs(cur_i).mMaxToxGrade= max_tox;
            CIobjs(cur_i).mFlgCensor = false; % complication w/in 6 m
            CIobjs(cur_i).mDateCensor = datenum(mFullToxDates(cur_inds(comp_ind))); % date of maximum toxicity
            CIobjs(cur_i).mDateComp = datenum(cur_max_tox_date); % complication w/in 6 m
            
            
        end
        
    end
    
    
    disp(['loading patient ',mPtIds{i}, ' was censored: ',num2str(CIobjs(cur_i).mFlgCensor)]);
    
    CIobjs(cur_i).mID = mPtIds{i};
    
    CIobjs(cur_i).mFxNum = cur_num_fractions;
    CIobjs(cur_i).mDoseTx = cur_tx;
    CIobjs(cur_i).mDateStartTx = datenum(cur_start_date(1));
    CIobjs(cur_i).mDateEndTx = datenum(cur_end_date(1));
    CIobjs(cur_i).mDateLastFollowUp = datenum(cur_followup_date(1));
    
    CIobjs(cur_i).mDoseBins_org = dosebins;
    %tmp
    %CIobjs(cur_i).mDoseBins_org = [dosebins(1);dosebins(2:5:end)];
    
    %# Get DVH info
    PtInfo.mLabel = cur_pt;
    PtInfo = PtInfo.fExtractRowData();
    flgdvh = PtInfo.mFlg;
    cvol = cell2mat(PtInfo.mData(flgdvh));
    
    cvol = cvol./cvol(1);
    CIobjs(cur_i).mVolCum = [cvol;0];
    CIobjs(cur_i) = CIobjs(cur_i).fCum2Diff();
    
    if ~strcmp(a2b_corr,'PHYS')
        CIobjs(cur_i).mBeta2AlphaCorrection = a2b_corr;
        CIobjs(cur_i).mBeta2Alpha = [1./alpha2beta];
    end
    cur_i=cur_i+1;
end

k=1;
CIobjs(cur_i:end)=[];
CGobjs(k) = CGobjs(k).fAddPatient(CIobjs);
CGobjs(k) = CGobjs(k).fCalculateEUD();


if to_save
    CGstrct = ObjToStruct(CGobjs);
    fn=['Z:/elw/MATLAB/meta_analy/meta_data/RTOG_',eud_binning,'_EUD_meta.mat'];
    save(fn,'Regions','CGobjs','CGstrct');
end
disp([]);
disp(['= Omitted patients =']);
disp(['  Missing DVHs: ', num2str(omitted_dvh)]);
disp(['  Excluded DVHs: ', num2str(excluded_dvh)]);
disp(['  Less than 6 months follow-up (w/o comp): ', num2str(omitted_lt6m)]);
disp(['  Etc: ', num2str(omitted_etc)]);
disp(['==='])
disp(['Total Number of patients: ',num2str(cur_i-1)]);%1 because cur_i starts at 1
disp(['Number of complications: ',num2str(n_comp)]);
toc;
end