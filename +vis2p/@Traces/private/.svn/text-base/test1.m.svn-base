
%% V2 FLD classifier with posterior

% variance normalization
vnorm  = @(x) bsxfun(@rdivide,x,std(x));

% get keys with V1-V2 recordings
keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

% intialize variables
[ad, dv, ed, rd, mn, ps,cpV2] = initialize('cell',length(keys),1);

parfor ikey = 1:length(keys);display(num2str(ikey)); key = keys(ikey);
    
    % set params
    key.trace_opt = 28;
    key.stats_opt =1;
    
    % get data
    tracesV1 = getTraces(StatArea(key),'area','v1o','key',key);
    tracesV2 = getTraces(StatArea(key),'area','v2o','key',key);
    if iscell(tracesV1)
        mbins = min(cellfun(@(x) size(x,2),tracesV1));
        mtrials = min(cellfun(@(x) size(x,3),tracesV1));
        tracesV1 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV1,'uniformoutput',0));
        tracesV2 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV2,'uniformoutput',0));
    end
    istart = 1000/fetch1(StatsSitesParams(key),'binsize')*1;
    iend = 1000/fetch1(StatsSitesParams(key),'binsize')*3;
    tracesV1 = tracesV1(:,istart:iend,:);
    tracesV2 = tracesV2(:,istart:iend,:);
    %     tracesV1 = tracesV2;
    mtraces = mean(tracesV1,3);
    
    % calculate posterior for V2 responses
    cp = nan(size(tracesV2,2),size(tracesV2,3));
    for iTrial = 1:size(tracesV2,3)
        % format data
        trIdx = 1:size(tracesV2,3);trIdx(iTrial) = [];
        train = tracesV2(:,:,trIdx);
        test = tracesV2(:,:,iTrial);
        
        % build the classifier
        group = reshape(bsxfun(@times,ones(size(train,2),size(train,3))',1:size(train,2))',[],1);
        obj = ClassificationDiscriminant.fit(train(:,:)',group);
        [~,score]= predict(obj,test(:,:)');
        cp(:,iTrial) =  diag(score);
    end
    cpV2{ikey} = cp(:);
    
    % compute V1 measures
    [rdist, edist,dvec, adist] = initialize('nan',size(tracesV2,2),size(tracesV2,3));
    for istim = 1:size(tracesV2,2)
        vcV1 = squeeze(tracesV1(:,istim,:));
        rdist(istim,:) = vnorm(pdist2(vcV1',mtraces(:,istim)','cosine'));
        edist(istim,:) = vnorm(pdist2(vcV1',mtraces(:,istim)'));
        dvec(istim,:) =  vnorm(dot(vcV1,repmat(mtraces(:,istim),[1 size(tracesV1,3)]))-norm(mtraces(:,istim)));
        adist(istim,:) = vnorm(pdist2(double(vcV1'),ones(1,size(mtraces,1)),'cosine'));
    end
    dv{ikey} = dvec(:);
    ed{ikey} = edist(:);
    rd{ikey} = rdist(:);
    ad{ikey} = adist(:);
    ps{ikey} = reshape(vnorm(squeeze(sparseness(tracesV1))')',[],1);
    mn{ikey} = reshape(vnorm(squeeze(sum(tracesV1))')',[],1);
end
save('28-1-data','ad', 'dv', 'ed', 'rd', 'mn', 'ps','cpV2')


%% classification and correlation structure 9:1

keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

[cr, crR,er, erR] = initialize('cell',length(keys),1);

parfor ikey = 1:length(keys);
    key = keys(ikey);
    
    key.trace_opt = 17;
    key.stats_opt =1;
    
    tracesV1 = getTraces(StatArea(key),'area','v1o','key',key);
    tracesV2 = getTraces(StatArea(key),'area','v2o','key',key);
    if iscell(tracesV1)
        mbins = min(cellfun(@(x) size(x,2),tracesV1));
        mtrials = min(cellfun(@(x) size(x,3),tracesV1));
        tracesV1 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV1,'uniformoutput',0));
        tracesV2 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV2,'uniformoutput',0));
    end
    istart = 1000/fetch1(StatsSitesParams(key),'binsize')*1;
    iend = 1000/fetch1(StatsSitesParams(key),'binsize')*3;
    tracesV1 = tracesV1(:,istart:iend,:);
    tracesV2 = tracesV2(:,istart:iend,:);
    
    % initialize
    [cr{ikey},er{ikey},crR{ikey},erR{ikey}] = initialize('nan',size(tracesV2,1),10);
    
    for igr = 1:10
        % format data
        trIdx = 1:size(tracesV1,3);tsIdx = trIdx(igr:10:end);trIdx(igr:10:end) = [];
        train = tracesV1(:,:,trIdx);
        test = tracesV1(:,:,tsIdx);
        trainGr = (tracesV2(:,:,trIdx)>0.005) + 1;
        testGr = (tracesV2(:,:,tsIdx)>0.005) + 1;
        
        
        % Iterate through all the cells in V2
        for iCell = 1:size(testGr,1)
            % equalize trials
            uniGr = unique(trainGr(iCell,:));
            [trainData, trainGroup, testData, testGroup,trainDataR] = ...
                initialize('cell',size(train,2),1);
            for iStim = 1:size(train,2)
                [tr, ts,trR,trG,tsG] = initialize('cell',length(uniGr),1);
                for iGr = 1:length(uniGr)
                    iUni = uniGr(iGr);
                    tr{iGr} = squeeze(train(:,iStim,squeeze(trainGr(iCell,iStim,:) == iUni)));
                    for iV1 = 1:size(tracesV1,1)
                        trR{iGr}(iV1,:) = tr{iGr}(iV1,randperm(size(tr{iGr},2)));
                    end
                    ts{iGr} = squeeze(test(:,iStim,squeeze(testGr(iCell,iStim,:) == iUni)));
                    trG{iGr} =  ones(1,size(tr{iGr},2))*iUni;
                    tsG{iGr} =  ones(1,size(ts{iGr},2))*iUni;
                end
                minTr = min(cellfun(@(x) size(x,2),tr));
                minTs = min(cellfun(@(x) size(x,2),ts));
                trainData{iStim}  = reshape(cell2mat(cellfun(@(x) x(:,1:minTr),...
                    tr,'UniformOutput',0)),size(tracesV1,1),[]);
                testData{iStim}  = reshape(cell2mat(cellfun(@(x) x(:,1:minTs),...
                    ts,'UniformOutput',0)),size(tracesV1,1),[]);
                trainDataR{iStim}  = reshape(cell2mat(cellfun(@(x) x(:,1:minTr),...
                    trR,'UniformOutput',0)),size(tracesV1,1),[]);
                trainGroup{iStim}  = reshape(cell2mat(cellfun(@(x) x(:,1:minTr),...
                    trG,'UniformOutput',0)),1,[]);
                testGroup{iStim}  = reshape(cell2mat(cellfun(@(x) x(:,1:minTs),...
                    tsG,'UniformOutput',0)),1,[]);
            end
            
            % with covariance matrix
            [C,err] = classify(cell2mat(testData')',cell2mat(trainData')',cell2mat(trainGroup')','quadratic');
            cr{ikey}(iCell,igr) = mean(C==cell2mat(testGroup')');
            er{ikey}(iCell,igr) = err;
            
            % without covariance matrix -randomize trials for each class defined from V2 response
            [C,err] = classify(cell2mat(testData')',cell2mat(trainDataR')',cell2mat(trainGroup')','quadratic');
            crR{ikey}(iCell,igr) = mean(C==cell2mat(testGroup')');
            erR{ikey}(iCell,igr) = err;
        end
    end
end


save('MeanClassifier-17-1','cr', 'er', 'ed', 'crR', 'erR')

%% classification and correlation structure 9:1

keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

[cr, crR,er, erR] = initialize('cell',length(keys),1);

parfor ikey = 1:length(keys);
    key = keys(ikey);
    
    key.trace_opt = 28;
    key.stats_opt =1;
    
    tracesV1 = getTraces(StatArea(key),'area','v1o','key',key);
    tracesV2 = getTraces(StatArea(key),'area','v2o','key',key);
    if iscell(tracesV1)
        mbins = min(cellfun(@(x) size(x,2),tracesV1));
        mtrials = min(cellfun(@(x) size(x,3),tracesV1));
        tracesV1 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV1,'uniformoutput',0));
        tracesV2 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV2,'uniformoutput',0));
    end
    istart = 1000/fetch1(StatsSitesParams(key),'binsize')*1;
    iend = 1000/fetch1(StatsSitesParams(key),'binsize')*3;
    tracesV1 = tracesV1(:,istart:iend,:);
    tracesV2 = tracesV2(:,istart:iend,:);
     trainGr = tracesV2>0.001;

    % initialize
    [cr{ikey},er{ikey},crR{ikey},erR{ikey}] = initialize('nan',size(tracesV2,1),size(tracesV2,2));
    for iStim = 1:size(tracesV1,2);
        for iCell = 1:size(tracesV2,1)
            
            testGr = squeeze(tracesV2(iCell,iStim,:)>0.001);
            if sum(testGr)<1 || sum(testGr)==size(tracesV2,3)
                continue
            end

            testi = {find(testGr)};
            testi{2} = find(~testGr);
             [~,imx]=max(cellfun(@length,testi));
              [~,imn]=min(cellfun(@length,testi));
            testi{imx} = testi{imx}(randperm(length(testi{imx}),length(testi{imn})));
            testGri = [testi{1};testi{2}];
            testGr = testGr(testGri) +1;
            test = squeeze(tracesV1(:,iStim,testGri));
            
            for igr = 1:length(testGri);
                tsI = testGri(igr);
                trI = 1:size(tracesV1,3);
                trI(tsI) = [];
               
                train = squeeze(tracesV1(:,iStim,trI));
                trainDataR = nan(size(train));
                for iV1 = 1:size(train,1)
                    trainDataR(iV1,:) = train(iV1,randperm(size(train,2)));
                end
                
                trainData  = train;
                testData  = test(:,igr);
                
                trainGroup  = squeeze(trainGr(iCell,iStim,trI));
                testGroup  = testGr(igr);
                
                try
                    
                    % with covariance matrix
                    [C,err] = classify(testData',trainData',trainGroup);
                    cr{ikey}(iCell,iStim,igr) = mean(C==testGroup');
                    er{ikey}(iCell,iStim,igr) = err;
                    
                    % without covariance matrix -randomize trials for each class defined from V2 response
                    [C,err] = classify(testData',trainDataR',trainGroup);
                    crR{ikey}(iCell,iStim,igr) = mean(C==testGroup');
                    erR{ikey}(iCell,iStim,igr) = err;
                end
            end
        end
    end
end

save('MeanClassifier-28-1','cr', 'er', 'ed', 'crR', 'erR')
%%
scatterhist1(1-  nanmean(cell2mat(er),2),1- nanmean(cell2mat(erR),2),...
    'names',{'FLD','FLD - w/o Covariance Structure'},...
    'title','Training performance - 500msec - PCA - 90 ')
scatterhist1( nanmean(cell2mat(cr),2),nanmean(cell2mat(crR),2),...
    'names',{'FLD','FLD - w/o Covariance Structure'},...
    'title','Testing Performance - 500msec -PCA -10')

%% v1 - V2 FLD classifier with posterior - noise correlations

% variance normalization
vnorm  = @(x) bsxfun(@rdivide,x,std(x));

% get keys with V1-V2 recordings
keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

% intialize variables
[cpV1,cprV1,cpV2,psV1,psrV1,psV2] = initialize('cell',length(keys),1);

parfor ikey = 1:length(keys);display(num2str(ikey)); key = keys(ikey);
    
    % set params
    key.trace_opt = 28;
    key.stats_opt =1;
    
    % get data
    tracesV1 = getTraces(StatArea(key),'area','v1o','key',key);
    tracesV2 = getTraces(StatArea(key),'area','v2o','key',key);
    if iscell(tracesV1)
        mbins = min(cellfun(@(x) size(x,2),tracesV1));
        mtrials = min(cellfun(@(x) size(x,3),tracesV1));
        tracesV1 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV1,'uniformoutput',0));
        tracesV2 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV2,'uniformoutput',0));
    end
    istart = 1000/fetch1(StatsSitesParams(key),'binsize')*1;
    iend = 1000/fetch1(StatsSitesParams(key),'binsize')*3;
    tracesV1 = tracesV1(:,istart:iend,:);
    tracesV2 = tracesV2(:,istart:iend,:);
    
    % calculate posterior for V2 responses
    cp = nan(size(tracesV2,2),size(tracesV2,3));ps =cp;
    for iTrial = 1:size(tracesV2,3)
        % format data
        trIdx = 1:size(tracesV2,3);trIdx(iTrial) = [];
        train = tracesV2(:,:,trIdx);
        test = tracesV2(:,:,iTrial);
        
        % build the classifier
        group = reshape(bsxfun(@times,ones(size(train,2),size(train,3))',1:size(train,2))',[],1);
        obj = ClassificationDiscriminant.fit(train(:,:)',group);
        [x,score]= predict(obj,test(:,:)');
        ps(:,iTrial) =  diag(score);
        cp(:,iTrial) = x==(1:size(test,2))';
    end
    psV2{ikey} = ps;
    cpV2{ikey} = cp;
    
    % calculate posterior for V1 responses
    cp = nan(size(tracesV1,2),size(tracesV1,3));ps =cp;
    for iTrial = 1:size(tracesV1,3)
        % format data
        trIdx = 1:size(tracesV1,3);trIdx(iTrial) = [];
        train = tracesV1(:,:,trIdx);
        test = tracesV1(:,:,iTrial);
        
        % build the classifier
        group = reshape(bsxfun(@times,ones(size(train,2),size(train,3))',1:size(train,2))',[],1);
        obj = ClassificationDiscriminant.fit(train(:,:)',group);
        [x,score]= predict(obj,test(:,:)');
        ps(:,iTrial) =  diag(score);
        cp(:,iTrial) = x==(1:size(test,2))';
    end
    psV1{ikey} = ps;
    cpV1{ikey} = cp;
    
    % calculate posterior for V1 responses without correlation structure
    cp = nan(size(tracesV1,2),size(tracesV1,3));ps =cp;
    for iTrial = 1:size(tracesV1,3)
        % format data
        trIdx = 1:size(tracesV1,3);trIdx(iTrial) = [];
        train = tracesV1(:,:,trIdx);
        test = tracesV1(:,:,iTrial);
        
        for iCell = 1:size(tracesV1,1);for iStim = 1:size(tracesV1,2); %#ok<*ALIGN>
            train(iCell,iStim,:) = train(iCell,iStim,randperm(size(train,3)));
        end;end
        
        % build the classifier
        group = reshape(bsxfun(@times,ones(size(train,2),size(train,3))',1:size(train,2))',[],1);
        obj = ClassificationDiscriminant.fit(train(:,:)',group);
        [x,score]= predict(obj,test(:,:)');
        
        ps(:,iTrial) =  diag(score);
        cp(:,iTrial) = x==(1:size(test,2))';
    end
    psrV1{ikey} = ps;
    cprV1{ikey} = cp;
    

end
save('28-1-NC V2V2','cpV1','cpV2','cprV1','psV1','psV2','psrV1')

%%
%% classification and correlation structure 9:1

keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

[cr, crR,er, erR] = initialize('cell',length(keys),1);

parfor ikey = 1:length(keys);
    key = keys(ikey);
    
    key.trace_opt = 17;
    key.stats_opt = 1;
    threshold = 0.001;
    
    tracesV1 = getTraces(StatArea(key),'area','v1o','key',key);
    tracesV2 = getTraces(StatArea(key),'area','v2o','key',key);
    if iscell(tracesV1)
        mbins = min(cellfun(@(x) size(x,2),tracesV1));
        mtrials = min(cellfun(@(x) size(x,3),tracesV1));
        tracesV1 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV1,'uniformoutput',0));
        tracesV2 = cell2mat(cellfun(@(x) x(:,1:mbins,1:mtrials),tracesV2,'uniformoutput',0));
    end
    istart = 1000/fetch1(StatsSitesParams(key),'binsize')*1;
    iend = 1000/fetch1(StatsSitesParams(key),'binsize')*3;
    tracesV1 = tracesV1(:,istart:iend,:);
    tracesV2 = tracesV2(:,istart:iend,:);
    trainGr = (tracesV2>threshold)+1;
    
    % initialize
    [cr{ikey},er{ikey},crR{ikey},erR{ikey}] = initialize('nan',size(tracesV2,1),size(tracesV2,2));
    for iStim = 1:size(tracesV1,2);
        for iCell = 1:size(tracesV2,1)
            
            % define groups
            testGr = squeeze(tracesV2(iCell,iStim,:)>threshold);
            
            % check for unbalanced conditions (groups heavily unequal)
            if sum(testGr)<1 || sum(testGr)==size(tracesV2,3);continue;end
            
            % take the minimum # trials randomly selected
            testi = {find(testGr)}; testi{2} = find(~testGr);
            [~,imx]=max(cellfun(@length,testi));
            [~,imn]=min(cellfun(@length,testi));
            testi{imx} = testi{imx}(randperm(length(testi{imx}),length(testi{imn})));
            
            % test data indices (fewer than real data with equal class #)
            testGri = [testi{1};testi{2}];
            
            % update the class from [0 1] to [1 2]
            testGr = testGr(testGri) +1;
            
            % test data (fewer than real data with equal class #)
            test = squeeze(tracesV1(:,iStim,testGri));
            
            % iterate though test trials (with equal number of classes)
            % (leave-one-out cross-validation)
            for igr = 1:length(testGri);
                
                % select the testing trial index
                tsI = testGri(igr);
                
                % select the training trials indices
                trI = testGri;
                trI(igr) = [];
                
                % select the training data
                train = squeeze(tracesV1(:,iStim,trI));
                
                % create the data with randomized trials
                trainDataR = nan(size(train));
                for iV1 = 1:size(train,1) % each V1 cell
                    % each group id
                        ttr = train(iV1,:)>threshold;
                        tdata = train(iV1,ttr);
                        trainDataR(iV1,ttr) = tdata(randperm(size(tdata,2)));
                        ttr = train(iV1,:)<=threshold;
                        tdata = train(iV1,ttr);
                        trainDataR(iV1,ttr) = tdata(randperm(size(tdata,2)));
                end
                
                trainData  = train;
                testData  = test(:,igr);
                
                trainGroup  = squeeze(trainGr(iCell,iStim,trI));
                testGroup  = testGr(igr);
                
                try
                    % with covariance matrix
                    [C,err] = classify(testData',trainData',trainGroup);
                    cr{ikey}(iCell,iStim,igr) = (C==testGroup');
                    er{ikey}(iCell,iStim,igr) = err;
                    
                    % without covariance matrix -randomize trials for each class defined from V2 response
                    [C,err] = classify(testData',trainDataR',trainGroup);
                    crR{ikey}(iCell,iStim,igr) = (C==testGroup');
                    erR{ikey}(iCell,iStim,igr) = err;
                end
            end
        end
    end
end

save('MeanClassifier-17-1','cr', 'er', 'ed', 'crR', 'erR')

