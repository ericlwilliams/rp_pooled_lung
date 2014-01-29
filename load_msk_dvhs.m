function load_msk_dvhs
%% Pt_MSK
% this function integrates all MSKCC data into one file (data object)
tic;
eud_binning = 'crs'; % med, crs, fine

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
    
    dvh_loc = 'Z:/elw/MATLAB/original_data/MSK/DVHs/whole/';
    %cd(dvh_loc);
    dvh_dir = dir(dvh_loc);
    dvh_dir = dvh_dir(3:end); %eliminate . and .. entries

    PtInfo = classDataFromXls();
    Regions = {'Whole'    'Ipsi'    'Contra'    'Superior'    'Inferior'    'Left'  'right'}';
    CGobjs = classOutcomeAnalysis();
    
    a2b_corr='NTD';
    a2b_val = 3;
    
    if isequal(eud_binning,'fine')
        CGobjs.mLymanN = 10.^(-1:0.01:1)';
    elseif isequal(eud_binning,'med')
        CGobjs.mLymanN = 10.^(-1:0.02:1)';
    elseif isequal(eud_binning,'crs')
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    else %default coarse
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    end

%% read xls file
    fn = 'Z:/elw/MATLAB/original_data/MSK/MSK_EUD.xls';
        % xlsData - array of structures with sheet name (e.g. 'Whole', 'Ispi',
    % 'Contra', etc.) and raw data
    xlsData=xlsFileRead(fn);

%% parse sheet (build CGobjs) for each region
   for k = 1:1
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
       
        %NumFx
         PtInfo.mLabel = 'NumFx';
        PtInfo = PtInfo.fExtractColData();
        flgcomp = PtInfo.mFlg;
        numfx = cell2mat(PtInfo.mData(flgcomp)); % censor info
         
        % Name
        PtInfo.mLabel = 'surname';
        PtInfo = PtInfo.fExtractColData();
        flgcomp = PtInfo.mFlg;
        ptname = [PtInfo.mData(flgcomp)]; % censor info
        ptname = strrep(ptname,',','');
        
          % Total dose
        PtInfo.mLabel = 'Tx';
        PtInfo = PtInfo.fExtractColData();
        flgtx = PtInfo.mFlg;
        tx = cell2mat(PtInfo.mData(flgtx));
        
        
      
        %% Get DVH information to calculate EUDs
        % Load cumulative and differential DVHs
        
        %% TODO
        %%%% [cur_dDVH,cur_cDVH]=import_ddvh_file(strcat(heart_data_loc,cur_name),patient_num);
        %%%  strcat(heart_data_loc,cur_name) ==== Z:\elw\MATLAB\original_data\MSK\heart\ARENA_HEART_ABSVOL.TXT

        if length(ptname)~= size(ptcode,1)
            disp(['Number of patients from data does not match number from DVHs']);
            return;
        end
        
        % transfer data to objects
        CIobjs = classOutcomeIndividual();
        CIobjs = repmat(CIobjs,size(ptcode)); % build matrix of cOIs 
        for n = 1:size(ptcode,1)
           
            % get current DVH
           dvh_name = dvh_dir(n).name; 
           cur_dvh = strcat(dvh_loc,dvh_name);
           [cur_dDVH,cur_cDVH,cur_dDVH_name,cur_cDVH_name]=import_ddvh_file(cur_dvh,n);
            
           %cvol = cur_cDVH(:,2);
           dvol = cur_dDVH(1:2:end,2);
           dvol = [dvol;0];
         
         
           %dosebins = (cur_cDVH(1:2:end-1,1)+cur_cDVH(2:2:end,1))./2;
           dosebins = cur_dDVH(1:2:end,1);
           dosebin_width = median(dosebins(2:end)-dosebins(1:end-1));
           
           dosebins = [dosebins;dosebins(end)+dosebin_width];
           dosebins = dosebins./100; %cGy -> Gy
           
           evalin('caller',['clear ',cur_dDVH_name]);
           evalin('caller',['clear ',cur_cDVH_name]);
                      
           CIobjs(n).mID=ptcode{n};
            CIobjs(n).mFlgCensor = ptcensor(n);
            CIobjs(n).mFxNum = numfx(n);
            CIobjs(n).mDoseBins_org = dosebins;        
            CIobjs(n).mVolDiff = dvol;
            CIobjs(n).mDoseTx = tx(n)/100;
            %CIobjs(n).mDoseBins_LQ = dosebins;
            CIobjs(n).mBeta2AlphaCorrection = a2b_corr;
            CIobjs(n).mBeta2Alpha = 1/a2b_val;
            
            
        end
        % save data
        CGobjs(k) = CGobjs(k).fAddPatient(CIobjs);
        CGobjs(k) = CGobjs(k).fCalculateEUD();

                
    end




% save result
    CGstrct = ObjToStruct(CGobjs);
    fn=['Z:/elw/MATLAB/meta_analy/meta_data/MSK_',eud_binning,'_EUD_meta_ppm.mat'];
    save(fn,'Regions','CGobjs','CGstrct');
toc;
end