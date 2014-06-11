% path = 'C:\Users\WORKHORSE\Google Drive\Paper\FabianAnalysis\data\';
load(['responses_ica_natural_GlobScale.mat'])

% set the mFR to ~ 1HZ
R = R(:,1:2000)*100;

% introduce noise to the V1 trial responses
Rt = poissrnd(repmat(R',[1 1 100]));

%% Simulate V2 traces
c = corr(R');
% create the V2 sampling indexes
iV2 = nan(101,size(Rt,1),size(Rt,3));
for i = 1:101
    [~,si] = sort(c(i,:),'descend');
    %     iV2(i,:,:) = squeeze(mean(Rt(:,randperm(size(R,1),4),:),2));
    iV2(i,:,:) = poissrnd(squeeze(mean(Rt(:,si(2:4),:),2)));
end
% iV2(iV2<std(iV2(:))/2) = 0;
% dR = bsxfun(@rdivide,iV2,sum(iV2>2*std(iV2(:)))+5);
dR = iV2;
dR(dR<std(dR(:))) = 0';
dR = dR.^3;

%% classification and correlation structure 9:1

tracesV1 = permute(Rt,[2 1 3]);
tracesV2 = dR;

% initialize
[cr,er,crR,erR] = initialize('cell',10,1);

parfor igr = 1:10

    display(num2str(igr))
    
    % format data
    trIdx = 1:size(tracesV1,3);tsIdx = trIdx(igr:10:end);trIdx(igr:10:end) = [];
    train = tracesV1(:,:,trIdx);
    test = tracesV1(:,:,tsIdx);
    trainGr = (tracesV2(:,:,trIdx)>0) + 1;
    testGr = (tracesV2(:,:,tsIdx)>0) + 1;
        
    % Iterate through all the cells in V2
    [cr{igr},er{igr},crR{igr},erR{igr}] = initialize('nan',size(testGr,1),1);
    for iCell = 1:size(testGr,1)
        try
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
            [C,err] = classify(cell2mat(testData')',cell2mat(trainData')',cell2mat(trainGroup')');
            cr{igr}(iCell) = mean(C==cell2mat(testGroup')');
            er{igr}(iCell) = err;
            
            % without covariance matrix -randomize trials for each class defined from V2 response
            [C,err] = classify(cell2mat(testData')',cell2mat(trainDataR')',cell2mat(trainGroup')');
            crR{igr}(iCell) = mean(C==cell2mat(testGroup')');
            erR{igr}(iCell) = err;
        end
    end
end
 save('SimulatedNoisecorrelations','cr','er','crR','erR')

%%
 scatterhist1(1-  mean(cell2mat(er),2),1- mean(cell2mat(erR),2),...
    'names',{'FLD','FLD - w/o Covariance Structure'},...
    'title','Training performance - 500msec - PCA - 90 ')
scatterhist1( mean(cell2mat(cr),2),mean(cell2mat(crR),2),...
    'names',{'FLD','FLD - w/o Covariance Structure'},...
    'title','Testing Performance - 500msec -PCA -10')
