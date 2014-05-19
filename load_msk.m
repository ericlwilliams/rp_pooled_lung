function load_msk
%% Pt_MSK
% this function integrates all MSKCC data into one file (data object)
% deprecated?? or for use with regions only, see load_msk_dvhs for pooled
% rp data
tic;

%% initialization
% PtInfo: 
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
    PtInfo = classDataFromXls();
    Regions = {'Whole'    'Ipsi'    'Contra'    'Superior'    'Inferior'    'Left'  'right'}';
    CGobjs = classOutcomeAnalysis();
    
    %coarse - ONLY coarse from data
    eud_binning = 'crs';
    
    if isequal(eud_binning,'fine')
        CGobjs.mLymanN = 10.^(-1:0.01:1)';
    elseif isequal(eud_binning,'med')
        CGobjs.mLymanN = 10.^(-1:0.02:1)';
    elseif isequal(eud_binning,'crs')
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    else %default coarse
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    end


    CGobjs = repmat(CGobjs,[7,1]);% one analysis set for each lung region

%% read xls file
    fn = 'Z:/elw/MATLAB/original_data/MSK/MSK_EUD.xls';
        % xlsData - array of structures with sheet name (e.g. 'Whole', 'Ispi',
    % 'Contra', etc.) and raw data
    xlsData=xlsFileRead(fn);

%% parse sheet (build CGobjs) for each region
    for k = 1:size(xlsData,1)
        PtInfo.mXlsRaw = xlsData(k).xlsRaw;% raw data for sheet k
        % MRN
        PtInfo.mLabel = 'MRN,'; % get MRN info from raw data stored in PtInfo
        PtInfo = PtInfo.fExtractColData();
        flgptcode = PtInfo.mFlg; % flags for location of data in column of PtInfo.mData
        ptcode = PtInfo.mData(flgptcode); % data under PtInfo.mLabel for sheet k
        f = cellfun('isclass', ptcode, 'char'); % transfer numeric pt code to string
        f=~f;
        ptcode(f) = cellfun(@(x) num2str(x), ptcode(f), 'UniformOutput', false);
        % Complication results
        PtInfo.mLabel = 'gt gd3';
        PtInfo = PtInfo.fExtractColData();
        flgcomp = PtInfo.mFlg;
        ptcensor = ~cell2mat(PtInfo.mData(flgcomp)); % censor info
       
      
        
        
        % EUD values
        PtInfo.mLabel = 'logn:';
            f = cellfun(@(x) strcmpi(x,PtInfo.mLabel), PtInfo.mXlsRaw);
            [m,n] = find(f); % location of the column
            if length(m)>1
                for kk = 2:length(m)
                    if ~isequal(PtInfo.mXlsRaw(:,n(1)),PtInfo.mXlsRaw(:,n(kk)))
                        warning(['There are ',num2str(length(m)),' cells containing the same column name "',PtInfo.mLabel,'", but the data are not consistant, pick the first column.']);
                        break;
                    end
                end
            elseif isempty(m)
                error(['Column name "',PtInfo.mLabel,'" not found, can not continue.']);
            end
        ptlgn = cell2mat(PtInfo.mXlsRaw(m(1)+1,n(1):n(1)+20))'; % cell array of logn values
        ptdose = cell2mat(PtInfo.mXlsRaw(m(1)+2:m(1)+79,n(1):n(1)+20)); % matrix of gEUDs for (pt,logn)

        % transfer data to objects
        CIobjs = classOutcomeIndividual();
        CIobjs = repmat(CIobjs,size(ptcode)); % build matrix of cOIs 
        for n = 1:size(ptcode,1)
            CIobjs(n).mID=ptcode{n};
            CIobjs(n).mFlgCensor = ptcensor(n);
            %CIobjs(n).mDoseTx = tx(n);
            
            CIobjs(n).mLymanN = ptlgn;
            CIobjs(n).mEUD = ptdose(n,:)';% gEUDs for each logn value
        end

        % save data
        CGobjs(k) = CGobjs(k).fAddPatient(CIobjs); % loads all (78) patient info, CGobjs basically list of CIobjs (classOutcomeIndividuals)
    end

%% Add left & right info
    fn = 'Z:/elw/MATLAB/regions/meta/Lt_rt_2004';
    [~,~,raw]=xlsread(fn,'Sheet1'); % read .xls sheet
    PtInfo.mXlsRaw = raw;
    % MRN
    PtInfo.mLabel = 'MRN';
    PtInfo = PtInfo.fExtractColData();
    flgptcode = PtInfo.mFlg;
    ptcode = PtInfo.mData(flgptcode); % cell array of MRN (PtInfo.mLabel) data
    f = cellfun('isclass', ptcode, 'char'); % transfer numeric pt code to string
    f=~f;
    ptcode(f) = cellfun(@(x) num2str(x), ptcode(f), 'UniformOutput', false);
    % GTV location
    PtInfo.mLabel = 'targetloc';
    PtInfo = PtInfo.fExtractColData();
    flgptcode = PtInfo.mFlg;
    loc = PtInfo.mData(flgptcode);
    % identifies left by looking for first character 'L' in column data
    flgl = cellfun(@(x) strcmpi(x(1),'L'),loc); 
    % assign left and right
    ptipsi = CGobjs(2).mGrp; % list of patients gEUD values for ipsi lung region
    ptcontra = CGobjs(3).mGrp;% list of patients gEUD values for contra lung region
    ptcodeipsi = {ptipsi.mID};
    ptcodecontra = {ptcontra.mID};
    if ~isequal(ptcodeipsi, ptcodecontra)
        error('patients in the ipsi and contra sheets are not in the same order');
    end
    ptl = ptcontra; ptr = ptipsi; % assume right lung is ipsi-lung
    f = ismember(ptcodeipsi,ptcode(flgl)); % find those whose ipsi-lung is left lung
    ptl(f) = ptipsi(f); ptr(f) = ptcontra(f); % make adjustment
    % save data
    CGobjs(k+1) = CGobjs(k+1).fAddPatient(ptl);% add new CGobjs for left pt info
    CGobjs(k+2) = CGobjs(k+2).fAddPatient(ptr);% add new CGobjs for right pt info
%% Add heart-target overlap information


% save result
    CGstrct = ObjToStruct(CGobjs);
    fn=['Z:/elw/MATLAB/meta_analy/meta_data/MSK_',eud_binning,'_EUD_meta_ppm.mat'];
    
    save(fn,'Regions','CGobjs','CGstrct');
toc;
end