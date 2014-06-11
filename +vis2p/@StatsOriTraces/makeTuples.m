function makeTuples( obj, key )

% get trace
trace = fetch1(Traces(key),'trace');
fps    = fetch1( Movies(key), 'fps' );

% Load times of the trace
time = fetch1(VisStims(key),'frame_timestamps')';

%stim times and types
stimTimes = fetchn(StatsPresents(key),'movie_times');
movieType = fetchn(StatsPresents(key),'movie_num');
movieTypes = bsxfun(@eq,movieType,unique(movieType)');

% get orientation stim info
R = fetch1(StatsOri(key),'R');
R = R(:,:,1:size(stimTimes{1}),1);

% find trace segments
traces = cell(1,length(stimTimes));
times = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = trace(time > stimTimes{iTimes}(1) & time < stimTimes{iTimes}(end));
    times{iTimes} = time(time > stimTimes{iTimes}(1) & time < stimTimes{iTimes}(end));
end

% equalize
traces = cellfun(@(x) x(1:min(cellfun(@length,traces))),traces,'UniformOutput',0);
times = cellfun(@(x) x(1:min(cellfun(@length,times))),times,'UniformOutput',0);

% collapse segments
traces = cell2mat(traces);
times = cell2mat(times);

% make sure there are no edge effects by removing 0.5 sec from the start and
% 0.5 sec from the end of the trials
traces = traces(round(fps/2):end - round(fps/2),:);
times = times(round(fps/2):end - round(fps/2),:);

%bin
d = max(1,round(key.binsize/1000*fps));
k = ones(d,1)/d;
trace = conv2(traces,k,'valid');
time = conv2(times,k,'valid');
if key.undersample
    trace = trace(1:d:end,:);
    time = time(1:d:end,:);
end

% bin the stim ori info
resp = nan(size(R,1),size(R,2),size(time,1));
for iTimes = 1:size(time,1)
    resp(:,:,iTimes) = mean(R(:,:,stimTimes{1} < time(iTimes,1) - key.resp_delay),3);
end

% equalize trials
extra = find(sum(movieTypes,1)>min(sum(movieTypes)));
indx = true(size(movieType,1),1);
for i = 1:length(extra)
    indx(find(movieTypes(:,extra(i)) == 1,1,'first')) = false;
end
movieType = reshape(movieType(indx),[min(sum(movieTypes)) length(unique(movieType))]);
trace = reshape(trace(:,indx)',min(sum(movieTypes)),length(unique(movieType)),[]);
[~, indx] = sort(movieType(1,:));
trace = trace(:,indx,:);

% do it
trace = reshape(trace,size(trace,1),[]);
resp = reshape(permute(resp,[2 1 3]),size(resp,2),[])';

ori = nan(size(trace,1),size(resp,2));
for iStim = 1:size(trace,1)
    ori(iStim,:) = corr(trace(iStim,:)',resp);
end

% make sure everything is greater than zero
baseline = min(ori(:));
ori = ori - baseline;

m = mean(ori);
v = var(ori);

% fit Von Mises with the true data
orientations = ((1:length(m))-1)*(pi/length(m));
FitVonMisses = fitVonMisesOne(m,orientations);
[~, Pdm] = opticalPropertiesONE(FitVonMisses);
key.Pdm = Pdm;

% Calciulate dprime
[~, dm] = max(m);
[dm dm90] = circOriOne(orientations,dm);
otiRaw = (m(dm)-m(dm90)) ./(m(dm)+m(dm90));
dPrimeOri = (m(dm)-m(dm90))/sqrt((v(dm)^2+v(dm90)^2)/2);

% Fit Von Mises for the shuffled data
binAreaMean = ori(:);
condIdxSfl = zeros(numel(binAreaMean),key.shuffle);
conditionIndexNum = numel(binAreaMean);


% shuffle the index
conditionIndexNumShuffled = 1:conditionIndexNum;
for i = 1:key.shuffle
    conditionIndexNumShuffled = conditionIndexNumShuffled(randperm(conditionIndexNum));
    condIdxSfl(:,i) = conditionIndexNumShuffled;
end

% generate randomly shuffled data
randBinArea = binAreaMean(condIdxSfl);
randBinArea = reshape(randBinArea,size(ori,1),size(ori,2),key.shuffle);

% mean the areas
areaMeanShfl = squeeze(mean(randBinArea,1));
varMeanShfl = squeeze(var(randBinArea,[],1));


% compute d' adn tuning indexes
[~, dm] = max(areaMeanShfl);
dm90 = zeros(1,key.shuffle);
for i = 1:key.shuffle
    [dm(i) dm90(i)] = circOriOne(orientations,dm(i));
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



