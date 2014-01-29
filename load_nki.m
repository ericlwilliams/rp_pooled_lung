function load_nki
% this function integrate all NKI data into one file
tic;
% prepare
    PtInfo = classDataFromXls();
    Regions = {'Whole'    'Ipsi'    'Contra'    'Superior'    'Inferior'    'Left'  'right'}';
    R_NKI = {'Whole'; 'ipsi'; 'contra'; 'cranial'; 'caudal'; 'Left'; 'Right'};
    CGobjs = classOutcomeAnalysis();
    
    eud_binning='crs';
        
    if isequal(eud_binning,'fine')
        CGobjs.mLymanN = 10.^(-1:0.01:1)';
    elseif isequal(eud_binning,'med')
        CGobjs.mLymanN = 10.^(-1:0.02:1)';
    elseif isequal(eud_binning,'crs')
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    else %default coarse
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    end

    
    CGobjs = repmat(CGobjs,[7,1]);

% read xls files
    if 1
        fn = 'Z:/elw/MATLAB/original_data/NKI/dvhregios.xls'; % regional DVS
        xlsData=xlsFileRead(fn);
        fn = 'Z:/elw/MATLAB/original_data/NKI/NSCLC_all_rev (2).xls'; % basic data
        xlsD = xlsFileRead(fn);
        save('Z:/elw/MATLAB/original_data/NKI/temp/tmp_nki.mat','xlsData','xlsD');
    else
        load('Z:/elw/MATLAB/original_data/NKI/temp/tmp_nki.mat','xlsData','xlsD');
    end

% parse basic data, choose pt with DE numbers
    f = cellfun(@(x) strcmp(x,'Basic data'), {xlsD.SheetName}); % find out the basic data sheet
    PtInfo.mXlsRaw = xlsD(f).xlsRaw;
    PtInfo.mLabel = 'DEnumber'; % find the pt with DE number
    PtInfo = PtInfo.fExtractColData();
    flgDE = PtInfo.mFlg; % pt with DE number are flaged
    flgDE(3) = false; % the third row is statistic result, not DE number
    PtInfo.mLabel = 'LFOnumber/Michigan number'; % find pt codes
    PtInfo = PtInfo.fExtractColData();
    ptcode = PtInfo.mData;
    ptcode = ptcode(flgDE); % pt id with DE number
    f = cellfun('isclass', ptcode, 'char');
    f=~f;
    ptcode(f) = cellfun(@(x) num2str(x), ptcode(f), 'UniformOutput', false);

    % RP info consistance check
    PtInfo.mLabel = 'Max RP grade'; % RP grade
    PtInfo = PtInfo.fExtractColData();
    ptrp = PtInfo.mData(flgDE);
    ptrp = cell2mat(ptrp);
    PtInfo.mLabel = 'RPgrade>=2'; % binary RP info
    PtInfo = PtInfo.fExtractColData();
    ptrp2 = PtInfo.mData(flgDE);
    ptrp2 = logical(cell2mat(ptrp2)); % 1--complication, 0--censored
    if ~isequal(ptrp>=2,ptrp2)
        error('patient complication info not consistent')
    end
    ptrp2 = ~ptrp2; % 1--censored, 0--complication

    
    
      % RP tx
    PtInfo.mLabel = 'Prescribed dose'; % RP grade
    PtInfo = PtInfo.fExtractColData();
    tx = PtInfo.mData(flgDE);
    tx = cell2mat(tx);
    
          % RP nfx
    PtInfo.mLabel = '#fractions (fraction size: 2.25 Gy'; % RP grade
    PtInfo = PtInfo.fExtractColData();
    nfx = PtInfo.mData(flgDE);
    nfx = cell2mat(nfx);
    
    
    
% parse whole lung DVHs
    CIobjs = classOutcomeIndividual();
    CIobjs = repmat(CIobjs,size(ptcode,1),1);

    f = cellfun(@(x) strcmp(x,'Whole'), {xlsD.SheetName}); % find out the whole lung sheet
    PtInfo.mXlsRaw = xlsD(f).xlsRaw;
    % pick up patients with DE numbers
    PtInfo.mLabel = 'LFOnumber';
    PtInfo = PtInfo.fExtractColData();
    f = cellfun('isclass',PtInfo.mData,'char'); % conver non-char element to char element to feed the ismember function
    f = ~f; % find non-character cells
    PtInfo.mData(f) = cellfun(@(x) num2str(x),PtInfo.mData(f),'UniformOutput', false);
    [flgpt, ptloc] = ismember(ptcode,PtInfo.mData);
    % DVH data location in the sheet
    PtInfo.mLabel = '#fractions (fraction size: 2.25 Gy';
    PtInfo = PtInfo.fExtractRowData();
    flgdvh = PtInfo.mFlg;
    dosebins = cell2mat(PtInfo.mData(flgdvh));
    %dosebins = [dosebins-0.5;0]; % -0.5 to shift the bin to the left boundary, add 0 to make sure the last bin corresponding to the zero volume
    
    dosebins = [dosebins;dosebins(end)+0.5]
    dosebins = dosebins-0.5; % -0.5 to shift the bin to the left boundary, add 0 to make sure the last bin corresponding to the zero volume
    
    % DVH for each patient
    f = find(flgpt); % patients with DE number
    flgpt(:)=false;
    for n = 1:length(f)
        vol = cell2mat(PtInfo.mXlsRaw(ptloc(f(n)),flgdvh))';
        if any(isnan(vol))
            continue;
        end
        CIobjs(n).mID = ptcode{f(n)}; % pt id
        CIobjs(n).mFlgCensor = ptrp2(f(n)); % censor info
        CIobjs(n).mDoseBins_org = dosebins; 
        CIobjs(n).mDoseBins_LQ = dosebins;
        CIobjs(n).mDoseTx = tx(n);
        CIobjs(n).mFxNum = nfx(n);
        CIobjs(n).mVolDiff = [vol;0];
        flgpt(n) = true;
    end
    % save data
    k=1;
    CGobjs(k) = CGobjs(k).fAddPatient(CIobjs(flgpt));
    CGobjs(k) = CGobjs(k).fCalculateEUD();
    
% parse regional DVHs
    for k = 2:size(R_NKI,1)
        f = cellfun(@(x) strcmp(x,R_NKI{k}), {xlsData.SheetName}); % find out the regional sheet
        PtInfo.mXlsRaw = xlsData(f).xlsRaw;
        % pick up patients with DE numbers
        PtInfo.mLabel = 'LFOnr';
        PtInfo = PtInfo.fExtractColData();
        f = cellfun('isclass',PtInfo.mData,'char'); % conver non-char element to char element to feed the ismember function
        f = ~f; % find non-character cells
        PtInfo.mData(f) = cellfun(@(x) num2str(x),PtInfo.mData(f),'UniformOutput', false);
        [flgpt, ptloc] = ismember(ptcode,PtInfo.mData);
        % DVH data location in the sheet
        PtInfo.mLabel = 'Dosebins(Gy):';
        PtInfo = PtInfo.fExtractRowData();
        flgdvh = PtInfo.mFlg;
        dosebins = cell2mat(PtInfo.mData(flgdvh));
        dosebins = [dosebins-0.5;0]; % -0.5 to shift the bin to the left boundary, add 0 to make sure the last bin corresponding to the zero volume
        % DVH for each patient
        f = find(flgpt); % patients with DE number
        flgpt(:)=false;
        for n = 1:length(f)
            vol = cell2mat(PtInfo.mXlsRaw(ptloc(f(n)),flgdvh))';
            if any(isnan(vol))
                continue;
            end
            CIobjs(n).mID = ptcode{f(n)}; % pt id
            CIobjs(n).mFlgCensor = ptrp2(f(n)); % censor info
            CIobjs(n).mDoseBins_org = dosebins;
            CIobjs(n).mDoseBins_LQ = dosebins;
            CIobjs(n).mVolDiff = [vol;0];
            flgpt(n) = true;
        end
        % save data
        CGobjs(k) = CGobjs(k).fAddPatient(CIobjs(flgpt));
        CGobjs(k) = CGobjs(k).fCalculateEUD();
    end

% save result

CGstrct = ObjToStruct(CGobjs);
fn=['Z:/elw/MATLAB/meta_analy/meta_data/NKI_',eud_binning,'_EUD_meta_ppm.mat'];
save(fn,'Regions','CGobjs','CGstrct');

toc;
end