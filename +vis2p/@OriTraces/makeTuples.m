function makeTuples( obj, key )


import vis2p.*

% get trace
trace = fetch1(Traces(key),'trace');
fps    = fetch1( Movies(key), 'fps' );

% Load times of the trace
[time, stim] = fetch1(VisStims(key),'frame_timestamps','stim_file');

%stim times and types
stimTimes = fetchn(OriPresents(key),'ori_times');
movieType = fetchn(OriPresents(key),'ori_num');
movieTypes = bsxfun(@eq,movieType,unique(movieType)');

% find trace segments
traces = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = trace(time' > stimTimes{iTimes}(1) & time' < stimTimes{iTimes}(end));
end

% get params 
[shuffle, resp_delay, resp_period] = fetch1(OriTracesParams(key),'shuffle','resp_delay','resp_period');

% equalize
indx = cellfun(@(x) length(x)>= ...
    floor(fps*mean(cellfun(@diff,stimTimes))/1000 - std(cellfun(@length,traces))),traces);
traces = traces(indx);
movieTypes = movieTypes(indx,:);
movieType = movieType(indx);
traces = cellfun(@(x) x(1:min(cellfun(@length,traces))),traces,'UniformOutput',0);

% collapse segments
trace = cell2mat(traces);

% equalize trials
extra = find(sum(movieTypes,1)>min(sum(movieTypes)));
indx = true(size(movieType,1),1);
for i = 1:length(extra)
    indx(find(movieTypes(:,extra(i)) == 1,sum(movieTypes(:,extra(i)),1) - ...
        min(sum(movieTypes)),'last')) = false;
end
movieType = reshape(movieType(indx),[min(sum(movieTypes)) length(unique(movieType))]);
trace = reshape(trace(:,indx)',min(sum(movieTypes)),length(unique(movieType)),[]);
[~, indx] = sort(movieType(1,:));
trace = trace(:,indx,:);

% do it
mx = min([size(trace,3) - round(resp_delay/1000*fps) round(resp_period/1000*fps)]);
ori = mean(trace(:,:,round(resp_delay/1000*fps):mx),3);

% make sure everything is greater than zero
baseline = min(ori(:));
ori = ori - baseline;

m = double(mean(ori,1));
v = double(var(ori,[],1));

% fit Von Mises with the true data
orientations = stim.params.constants.orientation/360 * 2 * pi;
FitVonMisses = fitVonMises(m,orientations);
[~, Pdm] = opticalProperties(FitVonMisses);
key.Pdm = Pdm;

% Calciulate dprime
[~, dm] = max(m);
[dm, ~,dm90] = circOri(orientations,dm);
otiRaw = (m(dm)-m(dm90)) ./(m(dm)+m(dm90));
dPrimeOri = (m(dm)-m(dm90))/sqrt((v(dm)^2+v(dm90)^2)/2);

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
varMeanShfl = squeeze(var(randBinArea,[],1));


% compute d' adn tuning indexes
[~, dm] = max(areaMeanShfl);
dm90 = zeros(1,shuffle);
for i = 1:shuffle
    [dm(i), ~, dm90(i)] = circOri(orientations,dm(i));
end

% convert to linear indicies
dm = sub2ind(size(areaMeanShfl),dm,1:size(areaMeanShfl,2));
dm90 = sub2ind(size(areaMeanShfl),dm90,1:size(areaMeanShfl,2));

pref = areaMeanShfl(dm);
orthPref = areaMeanShfl(dm90);
prefMorthPref = pref - orthPref;
varPref = varMeanShfl(dm);
oti = prefMorthPref ./(pref + orthPref);
dPrOriShfl = prefMorthPref ./...
    sqrt((varPref.^2 + (varMeanShfl(dm90)).^2)/2);

% Compute the significance
key.Poti = mean( oti > otiRaw );
key.Pdoti = mean( dPrOriShfl > dPrimeOri );
key.ori = ori + baseline;

% insert data
insert( obj, key );



