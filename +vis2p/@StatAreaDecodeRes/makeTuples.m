function makeTuples( obj, key )

% Get params
[area, class,trials,source,source_type] = fetch1(StatAreaDecodeResParams(key),...
    'area','area_classifier','trials','source','source_type');
k = [];
[k.ncells,k.nclasses,k.repetitions,k.classifier] = fetch1(StatAreaDecodeParams(key),...
    'ncells','nclasses','repetitions','classifier');

try
    % get source data
    if ~strcmp(source,'Stim')
        ks = key;
        ks.dec_opt = fetch1(StatAreaDecodeParams(k,['type="' source_type '"']),'dec_opt');
        cpV = fetch1(StatAreaDecode(ks),['cp' source]);
    end

    % Get area data
    k.type = fetch1(StatAreaDecodeParams(key),'type');
    k.classifier = class;
    ks = key;
    ks.dec_opt = fetch1(StatAreaDecodeParams(k),'dec_opt');
    cp = fetch1(StatAreaDecode(ks),['cp' area]);
catch
    display('Source not computed!')
    return
end

% Initialize
CMTOT = zeros(2,2);
nclasses = size(cp,3);
ntrials = size(cp,2);
npairs = size(cp,1);
mi = cell(size(npairs,1),1);

% Run for each pair
parfor ipair = 1:npairs
    
    % initialize
    F = zeros(nclasses);
    [CA,CR,FP,FN] = initialize('zeros',nclasses,1);
    
    % loop through trials
    for iTrial = 1:ntrials
        
        % loop through classes
        for iClass = 1:nclasses
            if strcmp(trials,'false');  indx1 = cpV(ipair,iTrial,iClass);
                if indx1==iClass; continue;end
            elseif strcmp(trials,'true');   indx1 = cpV(ipair,iTrial,iClass);
                if indx1~=iClass; continue;end
            end
            
            if ~strcmp(source,'Stim')
                i = cpV(ipair,iTrial,iClass);
            else
                i = iClass;
            end
            indx = cp(ipair,iTrial,iClass);
            F(i,indx) = F(i,indx) + 1;
        end
    end
    
    % Calculate confusion matrix
    d = diag(F,0);
    for iclass = 1:nclasses
        CA(iclass) = F(iclass,iclass);
        dind = true(size(d));dind(iclass) = false;
        CR(iclass) = sum(d(dind));
        FN(iclass) = sum(F(iclass,dind));
        FP(iclass) = sum(F(dind,iclass));
    end
    CM = zeros(2,2);
    CM(1,1) = sum(CA);
    CM(1,2) = sum(FN);
    CM(2,1) = sum(FP);
    CM(2,2) = sum(CR);
    
    % update confusion matrix
    CMTOT = CMTOT + CM;
    
    % calculate mutual information
    pA = CM/sum(CM(:));
    pi = sum(CM,2)/sum(CM(:));
    pj = sum(CM,1)/sum(CM(:));
    pij = pi*pj;
    if (sum(FN)+sum(FP)) == 0 && (sum(CA)+sum(CR))>0 % this is wrong, it should be FN+FP
        mi{ipair} = 1;
    elseif (sum(CA)+sum(CR)) == 0 && (sum(FN)+sum(FP))>0
        mi{ipair} = 0;
    else
        mi{ipair} = sum(sum(pA.*log2(pA./pij)));
    end
end
miSite = nanmean(cell2mat(mi));

% calculate classification performance
if ~strcmp(source,'Stim')
    perf = sum(reshape(cp==cpV,[],1))/numel(cp);
else
    perf = sum(reshape(cp==repmat(permute(1:size(cp,3),...
        [1 3 2]),[size(cp,1),size(cp,2),1]),[],1))/numel(cp);
end

% save the data
key.cperf = perf;
key.mutinfo = miSite;
key.confmat = CMTOT;

% insert data
insert( obj, key );
