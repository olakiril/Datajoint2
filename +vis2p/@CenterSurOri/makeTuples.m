function makeTuples( obj, key )

import vis2p.*

% params
shuffle = 5000;
names = {'Out','In'};

% get traces [oriOut oriIn cells trials]
[traces,oris] = getTraces(CenterSurOri,'key',key,'compute',1);

% control error
if isempty(traces)
    display('Nothing to compute!')
    return
end

% Calculate tuning 
for iArea = 1:2
    idx = oris{iArea}>=0;
    if iArea==1
        ori = squeeze(traces(idx,~idx,:,:))';
    else
        ori= squeeze(traces(~idx,idx,:,:))';
    end
    
    % make sure everything is greater than zero
    baseline = min(ori(:));
    ori = ori - baseline;
    
    m = double(mean(ori,1));
    
    % fit Von Mises with the true data
    orientations = oris{iArea}(idx)'/360 * 2 * pi;
    FitVonMisses = fitVonMises(m,orientations);
    [~,~,~,Pdm] = opticalProperties(FitVonMisses);
    
    % Calciulate 
     [~, dm] = max(m);
%     [~, dm] = min(abs(Pdm*360 / 2 / pi - oris{iArea}(idx)));
    [dm,~,dm90,dm270] = circOri(orientations,dm);
    otiRaw = (m(dm)-mean(m([dm90 dm270]))) ./(m(dm)+mean(m([dm90 dm270])));
        
    % Fit Von Mises for the shuffled data
    binAreaMean = ori(:);
    condIdxSfl = zeros(numel(binAreaMean),shuffle);
    conditionIndexNum = numel(binAreaMean);
    
    % shuffle the index
    conditionIndexNumShuffled = 1:conditionIndexNum;
    for i = 1:shuffle
        conditionIndexNumShuffled = conditionIndexNumShuffled(randperm(conditionIndexNum));
        condIdxSfl(:,i) = conditionIndexNumShuffled;
    end
    
    % generate randomly shuffled data
    randBinArea = binAreaMean(condIdxSfl);
    randBinArea = reshape(randBinArea,size(ori,1),size(ori,2),shuffle);
    
    % mean the areas
    areaMeanShfl = squeeze(mean(randBinArea,1));
    
    % compute tuning indexes
    [~, dm] = max(areaMeanShfl);
    dm90 = zeros(1,shuffle);dm270 = dm90;
    for i = 1:shuffle
        [dm(i), ~,dm90(i),dm270(i)] = circOri(orientations,dm(i));
    end
    
    % convert to linear indicies
    dm = sub2ind(size(areaMeanShfl),dm,1:size(areaMeanShfl,2));
    dm90 = sub2ind(size(areaMeanShfl),dm90,1:size(areaMeanShfl,2));
    dm270 = sub2ind(size(areaMeanShfl),dm270,1:size(areaMeanShfl,2));
    
    pref = areaMeanShfl(dm);
    orthPref = mean([areaMeanShfl(dm90); areaMeanShfl(dm270)]);
    prefMorthPref = pref - orthPref;
    oti = prefMorthPref ./(pref + orthPref);

    % Compute the significance
    eval(['key.Poti' names{iArea} ' = mean( oti > otiRaw );'])
    eval(['key.ori' names{iArea} ' = ori + baseline;'])
    eval(['key.Pdm' names{iArea} ' = Pdm;'])
    
    if iArea==2
       key.fitVM = FitVonMisses; 
    end
end

% insert data
insert( obj, key );