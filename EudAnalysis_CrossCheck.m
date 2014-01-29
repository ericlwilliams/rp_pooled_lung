function EudAnalysis_CrossCheck
tic; %close all;


fig_loc = 'Z:/elw/MATLAB/meta_analy/figures/latest/';
%fn = 'Z:/elw/MATLAB/meta_analy/meta_data/MSK_NKI_UMich_RTOG_EUD_meta.mat';
fn = 'C:\Documents and Settings\williae1\meta_data\meta_lung\MSK_NKI_UMich_RTOG_EUD_meta.mat';

load(fn,'CGmsk','CGnki','CGum','CGrtog','CGcomb');

% Want bins of log10(a), but calc done with a
loga = -1:0.1:1;
a = 10.^loga;

CGs = {CGnki CGum CGrtog};
colors = {'b' 'm' 'g'};
figure;
for i=1:length(CGs)
    
    OCobj=CGs{i};
    grp = OCobj.mGrp;
    n_grp = OCobj.mNumInGrp;
    
    
    
    dosebins = grp(1).mDoseBins_LQ; % # vol bins X 1
    vols = [grp.mVolDiff]'; % # pts X # vol bins 
    if isempty(vols)
        for m=1:n_grp
            grp(m) = grp(m).fCum2Diff(); 
        end
       vols = [grp.mVolDiff]'; % # pts X # vol bins 
    end
    
    euds = zeros(n_grp,length(loga));
    for j=1:length(a)
        cur_a = a(j);
        for k=1:n_grp % loop of patients
            cur_eud = vols(k,:)*(dosebins.^(cur_a));
            cur_eud = cur_eud^(1/cur_a);
            euds(k,j) = cur_eud;
        end
    end
    
    % fit log reg
    pttotal = ones(n_grp,1);
    ptcomp = ~[grp.mFlgCensor]';
    pvals = ones(length(loga),1);
    llhds = ones(length(loga),1);
    for l=1:length(loga)
        [~,dev,s]=glmfit(euds(:,l),[ptcomp pttotal],'binomial','link','logit');
        pvalue = [s.p];
        pvals(l) = pvalue(2,:);
        
        llhds(l) = -0.5*(dev/s.dfe); % deviations
     
    end
    
    figure(1);hold on;
    semilogy(loga,pvals,'Color',colors{i});
    cur_ylim = [ylim];
    ylim([min(cur_ylim(1),min(pvals)) 1]);

    
    figure(2);hold on;
    plot(loga,llhds,'Color',colors{i});
    
    
end


end

