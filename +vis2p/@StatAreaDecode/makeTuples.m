function makeTuples( obj, key )

tuple = key;
tuple = rmfield(tuple,'stim_idx');
tuple = rmfield(tuple,'trial_trigger');
discrimType = 'linear';

% get traces
trV1 = getTraces(StatArea(key),'area','v1o','key',key,'compute2',1,'collapse',1);
trV2 = getTraces(StatArea(key),'area','v2o','key',key,'compute2',1,'collapse',1);

% control error
if isempty(trV1)
    display('Nothing to compute!')
    return
end

% get info
[nclasses, classifier, type] = fetch1(StatAreaDecodeParams(key),'nclasses','classifier','type');

switch type
    case 'rand' % randomize traces
        rindx = 1:numel(trV1);
        for i = 1:10000; rindx = rindx(randperm(numel(trV1))); end
        trV1 = reshape(trV1(rindx),size(trV1,1),size(trV1,2),size(trV1,3));
        rindx = 1:numel(trV2);
        for i = 1:10000; rindx = rindx(randperm(numel(trV2))); end
        trV2 = reshape(trV2(rindx),size(trV2,1),size(trV2,2),size(trV2,3));
        
    case 'cov' % get rid of covariance structure by randomizing trials
        trV1 = randdim(trV1,3);
        trV2 = randdim(trV2,3);
        discrimType = 'diagLinear';
        
    case 'var' % randomize variance across mean
        rindx = 1:numel(trV1);
        for i = 1:1000; rindx = rindx(randperm(numel(trV1))); end
        pm = repmat(mean(trV1,3),[1 1 size(trV1,3)]);
        vr = trV1(:) - pm(:);
        trV1 = pm + reshape(vr(rindx),size(trV1,1),size(trV1,2),size(trV1,3));
        rindx = 1:numel(trV2);
        for i = 1:1000; rindx = rindx(randperm(numel(trV2))); end
        pm = repmat(mean(trV2,3),[1 1 size(trV2,3)]);
        vr = trV2(:) - pm(:);
        trV2 = pm + reshape(vr(rindx),size(trV2,1),size(trV2,2),size(trV2,3));
end

% get the sizes
ntrials = size(trV1,3);
if nclasses==0;nclasses = size(trV1,2);end
pairs = nchoosek(1:size(trV1,2),nclasses);

% initialize
% cpV1 = ones(size(pairs,1),size(trV1,3),nclasses,'int8');
% cpV2 = ones(size(pairs,1),size(trV1,3),nclasses,'int8');
cpV1 = cell(size(pairs,1),1);
cpV2 = cell(size(pairs,1),1);

% loop through the pairs
parfor ipair = 1:size(pairs,1)
    %     try
    dataV1 = trV1(:,pairs(ipair,:),:);
    dataV2 = trV2(:,pairs(ipair,:),:);
    
    % loop through trials
    for iTrial = 1:ntrials
        
        % calculate mean without taking that trial into account
        ind = true(ntrials,1);
        ind(iTrial) = false;
        mV1 = nanmean(dataV1(:,:,ind),3);
        mV2 = nanmean(dataV2(:,:,ind),3);
        
        if strcmp(classifier,'nnclassRaw')
            % loop through classes
            for iClass = 1:nclasses
                [~,cpV1{ipair}(iTrial,iClass)] = min(pdist2(mV1',dataV1(:,iClass,iTrial)'));
                [~,cpV2{ipair}(iTrial,iClass)] = min(pdist2(mV2',dataV2(:,iClass,iTrial)'));
            end
        elseif strcmp(classifier,'SVM')
            mV1 = (dataV1(:,:,ind));
            mV2 = (dataV2(:,:,ind));
            groupV1 = reshape(bsxfun(@times,ones(size(mV1,2),size(mV1,3))',1:size(mV1,2))',[],1);
            groupV2 = reshape(bsxfun(@times,ones(size(mV2,2),size(mV2,3))',1:size(mV2,2))',[],1);
            
            SVMV1 = svmtrain(mV1(:,:)',groupV1(:));
            SVMV2 = svmtrain(mV2(:,:)',groupV2(:));
            
            % loop through classes
            for iClass = 1:nclasses
                cpV1{ipair}(iTrial,iClass) = svmclassify(SVMV1,dataV1(:,iClass,iTrial)');
                cpV2{ipair}(iTrial,iClass) = svmclassify(SVMV2,dataV2(:,iClass,iTrial)');
            end
        elseif strcmp(classifier,'fisher')
            
            mV1 = (dataV1(:,:,ind));
            mV2 = (dataV2(:,:,ind));
            groupV1 = reshape(bsxfun(@times,ones(size(mV1,2),size(mV1,3))',1:size(mV1,2))',[],1);
            groupV2 = reshape(bsxfun(@times,ones(size(mV2,2),size(mV2,3))',1:size(mV2,2))',[],1);
            try
                fldV1 = ClassificationDiscriminant.fit(mV1(:,:)',groupV1(:),'discrimType',discrimType);
                fldV2 = ClassificationDiscriminant.fit(mV2(:,:)',groupV2(:),'discrimType',discrimType);
            catch
                discrimType2 = 'pseudoLinear';
                fldV1 = ClassificationDiscriminant.fit(mV1(:,:)',groupV1(:),'discrimType',discrimType2);
                fldV2 = ClassificationDiscriminant.fit(mV2(:,:)',groupV2(:),'discrimType',discrimType2);
            end
            
            % loop through classes
            for iClass = 1:nclasses
                cpV1{ipair}(iTrial,iClass) = predict(fldV1,dataV1(:,iClass,iTrial)');
                cpV2{ipair}(iTrial,iClass) = predict(fldV2,dataV2(:,iClass,iTrial)');
            end
        end
    end
    %     catch err
    %        cpV1{ipair} = err.identifier;
    %     end
end


% reliability
tuple.cpV1 = reshape(cat(2,cpV1{:}),[],size(cpV1{1},1),size(cpV1{1},2));
tuple.cpV2 = reshape(cat(2,cpV2{:}),[],size(cpV2{1},1),size(cpV2{1},2));

% insert data
insert( obj, tuple );

