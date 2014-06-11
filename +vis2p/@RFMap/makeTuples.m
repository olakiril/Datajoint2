function makeTuples( obj, key )

tuple = key;

import vis2p.*

% Load trace
trace = fetch1(Traces(key),'trace');

% Load times of the trace
times = fetch1(VisStims(key),'frame_timestamps');

% Load Stimulus info
[dotTimes, dotLocX, dotLocY, dotColors] = ...
    fetchn(RFPresents(key),'dot_time','dot_loc_x','dot_loc_y','dot_color');

% get the specific time window
ontime = fetch1(RFOpts(key),'time_on');
offtime = fetch1(RFOpts(key),'time_off');
fractData = fetch1(RFOpts(key),'fraction_of_data');

% get the unique stimuli
uniX = unique(dotLocX);
uniY = unique(dotLocY);
uniColor = unique(dotColors);

% initialize bined trace
bTrace = nan(length(dotTimes),1);
X = nan(length(dotTimes),1);
Y = nan(length(dotTimes),1);
Color = nan(length(dotTimes),1);
keep = false(length(dotTimes),1);

% loop through all stimuli and find mean activity
for iStims = 1:length(dotTimes)

    % Create Time Windows
    binStart = dotTimes(iStims) + ontime;
    binEnd = dotTimes(iStims) + offtime;

    if  times(end)>= binEnd
        % Find the times for the windows
        binStartFrames = find(times>=binStart,1,'first');
        binEndFrames = find(times>=binEnd,1,'first');

        % Find the area
        bTrace(iStims) = mean(trace(binStartFrames:binEndFrames));
        keep(iStims) = true;
    end

    % find stim indx
    X(iStims) = find(uniX == dotLocX(iStims));
    Y(iStims) = find(uniY == dotLocY(iStims));
    Color(iStims) = find(uniColor == dotColors(iStims));
end

% take only valid samples
bTrace = bTrace(keep);

% find unique combinations.reverse axis for matrix indexing
allStims = [Y(keep) X(keep) Color(keep)];

% get the unique stimuli
uniX = unique(X(keep));
uniY = unique(Y(keep));
uniColor = unique(Color(keep));
uni = unique(allStims,'rows');

% initialize response matrix. reverse axis for matrix indexing
response = nan(length(uniY),length(uniX),length(uniColor),length(dotTimes));

% fill response matrix
for iUni = 1:length(uni)
    indx = find(sum(allStims == repmat(uni(iUni,:),[size(allStims,1) 1]),2)==3);
    response(uni(iUni,1),uni(iUni,2),uni(iUni,3),1:length(indx)) = bTrace(indx);
end

% % keep only full repetitions
% nonnan = find(isnan(sum(sum(sum(response,1),2),3)),1,'first') - 1;
% response = response(:,:,:,1:nonnan);

% get unique random seed and randomize data
seed = key2num(key);
rand('state',seed);
randIndx = randperm(size(response,4));

% select number of responses ( no less than 1 )
tuple.repeat_avg = max(1,round(size(response,4)* fractData));
dataIndx = randIndx(1:tuple.repeat_avg);

if length(uniColor)==1
    response(:,:,2,:,:) = response(:,:,1,:,:);
end

% get the trial averages
tuple.off_rf = nanmedian(response(:,:,1,dataIndx),4);
tuple.off_rf(isnan(tuple.off_rf)) = 0;
tuple.on_rf = nanmedian(response(:,:,2,dataIndx),4); % reverse for color indication 255 is on
tuple.on_rf(isnan(tuple.on_rf)) = 0;
tuple.onmoff_rf = nanmedian(response(:,:,2,dataIndx),4) - nanmedian(response(:,:,1,dataIndx),4);
tuple.onmoff_rf(isnan(tuple.onmoff_rf)) = 0;
tuple.onpoff_rf = nanmedian(response(:,:,2,dataIndx),4) + nanmedian(response(:,:,1,dataIndx),4);
tuple.onpoff_rf(isnan(tuple.onpoff_rf)) = 0;

insert( obj, tuple );






