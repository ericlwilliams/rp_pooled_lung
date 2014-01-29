function load_umich
% this function integrate all NKI data into one file
tic;
% prepare
    PtInfo = classDataFromXls();
    CGobjs = classOutcomeAnalysis();
    
    eud_binning='fine';
        
    if isequal(eud_binning,'fine')
        CGobjs.mLymanN = 10.^(-1:0.01:1)';
    elseif isequal(eud_binning,'med')
        CGobjs.mLymanN = 10.^(-1:0.02:1)';
    elseif isequal(eud_binning,'crs')
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    else %default coarse
        CGobjs.mLymanN = 10.^(-1:0.1:1)';
    end

    
    
    alpha2beta = 3;
    a2b_corr = 'NTD';
    CGobjs.mBeta2Alpha = [1./alpha2beta]; 
% read xls files
    fn = 'Z:/elw/MATLAB/original_data/UMich/MichLungData9204_Final.xls'; % basic data
    xlsD = xlsFileRead(fn);

% parse basic data, choose pt with DE numbers
    f = cellfun(@(x) strcmp(x,'Data'), {xlsD.SheetName}); % find out the basic data sheet
    PtInfo.mXlsRaw = xlsD(f).xlsRaw;
    PtInfo.mLabel = 'ID#'; % find the pt with DE number
    PtInfo = PtInfo.fExtractColData();
    flgDE = PtInfo.mFlg; % pt with DE number are flaged
    ptcode = PtInfo.mData;
    ptcode = ptcode(flgDE); % pt id with DE number
    f = cellfun('isclass', ptcode, 'char');
    f=~f;
    ptcode(f) = cellfun(@(x) num2str(x), ptcode(f), 'UniformOutput', false);

    % RP info 
    PtInfo.mLabel = 'Lung Toxicity Grade (SWOG w/ clarifications - see paper)'; % RP grade
    PtInfo = PtInfo.fExtractColData();
    ptrp = PtInfo.mData(flgDE);
    ptrp = cell2mat(ptrp);
    
    % RTOG RP >= Grade 2
    ptrp2 = ptrp>=2;
    ptrp2 = ~ptrp2; % 1--censored, 0--complication
    
    % Num FX
    PtInfo.mLabel = '# Fx';
    PtInfo = PtInfo.fExtractColData();
    numfx = PtInfo.mData(flgDE);
    numfx = cell2mat(numfx);
    
     % total dose
    PtInfo.mLabel = 'Total Dose [Gy]';
    PtInfo = PtInfo.fExtractColData();
    tx = PtInfo.mData(flgDE);
    tx = cell2mat(tx);
    
     % Days to fu
    PtInfo.mLabel = 'Elapsed days between end of RT and last FU';
    PtInfo = PtInfo.fExtractColData();
    datefu = PtInfo.mData(flgDE);
    datefu = cell2mat(datefu);
    
    % Days to comp
    PtInfo.mLabel = 'Elapsed days between RT start and G2 ';
    PtInfo = PtInfo.fExtractColData();
    datestart2comp = PtInfo.mData(flgDE);
    datestart2comp = cell2mat(datestart2comp);
    datestart2comp(isnan(datestart2comp))=Inf;
    
    PtInfo.mLabel = 'Elapsed days between RT start and end';
    PtInfo = PtInfo.fExtractColData();
    datestart2end = PtInfo.mData(flgDE);
    datestart2end = cell2mat(datestart2end);
    
    % date to complication from end of treatment
    date2comp = datestart2comp - datestart2end;
    
    
    
    
% parse whole lung DVHs
    CIobjs = classOutcomeIndividual();    
    CIobjs = repmat(CIobjs,size(ptcode,1),1);

    % DVH data location in the sheet
    PtInfo.mLabel = 'dosebins';
    PtInfo = PtInfo.fExtractRowData();
    flgdvh = PtInfo.mFlg;
    dosebins = cell2mat(PtInfo.mData(flgdvh));
    dosebins = [dosebins;dosebins(end)+0.5];
    dosebins = [dosebins-0.5]; % -0.5 to shift the bin to the left boundary, add 0 to make sure the last bin corresponding to the zero volume
    % DVH for each patient
    f = find(flgDE); % patients with DE number
    flgpt=flgDE;
    flgpt(:)=false;
    omitted_pts=0;
    tot_pts=0;
    comp_pts=0;

    for n = 1:length(f)
        
        % no complication if comp after 6 months
        if ~ptrp2(n) && date2comp(n)>180 % not acute
            ptrp2(n)=1; % don't count as complication
        end
            
        % requre at least 6 months of followup
        if datefu(n)<=180 && ptrp2(n) % followup less than 6 months and no comp
            omitted_pts=omitted_pts+1;
            disp(['Omitting pt ',num2str(ptcode{n}),10,...
                ' Date FU: ',num2str(datefu(n)),10,...
                ' Date to RP: ',num2str(date2comp(n))]);
            continue;
            
        end
        
       
                
        vol = cell2mat(PtInfo.mXlsRaw(f(n),flgdvh))';
        vol = vol./sum(vol);% to relative DVH
        if any(isnan(vol))
            omitted_pts=omitted_pts+1;
            disp(['Omitting pt ',num2str(ptcode{n}),10,...
                ' No DVH']);
            continue;
        end
        CIobjs(n).mID = ptcode{n}; % pt id
        CIobjs(n).mFlgCensor = ptrp2(n); % censor info
        CIobjs(n).mDoseBins_org = dosebins; 
        CIobjs(n).mDoseBins_LQ = dosebins;
        CIobjs(n).mFxNum = numfx(n);

        CIobjs(n).mVolDiff = [vol;0];
        
        CIobjs(n).mDoseTx = tx(n);
        flgpt(n) = true;
        
        % UMich data has LQ correction with alpha2beta = 2.5
        % undo this correction
        CIobjs(n) = CIobjs(n).fUndoLinearQuadratic(2.5);
        
        % redo LQ correction with alpha2beta=3
        CIobjs(n).mBeta2AlphaCorrection = a2b_corr;
        CIobjs(n).mBeta2Alpha = [1./alpha2beta];

        
        
        
        comp_pts=comp_pts+~ptrp2(n);
        tot_pts=tot_pts+1;
       
    end
    % save data
    disp([10,10,'Total patients: ',num2str(tot_pts),10,...
        'Omitted: ',num2str(omitted_pts),' w/ Comp: ',num2str(comp_pts)]);
    CGobjs = CGobjs.fAddPatient(CIobjs(flgpt));
    CGobjs = CGobjs.fCalculateEUD();

    
% save result

CGstrct = ObjToStruct(CGobjs);
fn=['Z:/elw/MATLAB/meta_analy/meta_data/UMich_',eud_binning,'_EUD_meta_ppm.mat'];
save(fn,'CGobjs','CGstrct');

toc;
end