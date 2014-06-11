function makeTuples( obj, key )

tuple = key;
rz = 0.5;

% Load trace
tp = tpReader(Scans(key));
trace =tp.imCh{1}(:,:,:);
trace = imresize(trace,rz);
sz = size(trace);
trace = permute(trace,[3 1 2]);
trace = trace(:,:);

% calc dfof
trace = bsxfun(@rdivide,double(trace),mean(trace)) - 1;

% Load times of the trace
times = fetch1(VisStims(key),'frame_timestamps');

% Load Stimulus info
[dotTimes, dotLocX, dotLocY, dotColors] = ...
    fetchn(RFPresents(key),'dot_time','dot_loc_x','dot_loc_y','dot_color');

% get the specific time window
ontime = fetch1(RFOpts(key),'time_on');
offtime = fetch1(RFOpts(key),'time_off');

% get the unique stimuli
uniX = unique(dotLocX);
uniY = unique(dotLocY);
uniColor = unique(dotColors);

% initialize bined trace
bTrace = nan(length(dotTimes),size(trace,2));
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
        bTrace(iStims,:) = mean(trace(binStartFrames:binEndFrames,:));
        keep(iStims) = true;
    end

    % find stim indx
    X(iStims) = find(uniX == dotLocX(iStims));
    Y(iStims) = find(uniY == dotLocY(iStims));
    Color(iStims) = find(uniColor == dotColors(iStims));
end

% take only valid samples
bTrace = bTrace(keep,:);

% find unique combinations.reverse axis for matrix indexing
allStims = [Y(keep) X(keep) Color(keep)];

% get the unique stimuli
uniX = unique(X(keep));
uniY = unique(Y(keep));
uniColor = unique(Color(keep));
uni = unique(allStims,'rows');

% initialize response matrix. reverse axis for matrix indexing
response = nan(length(uniY),length(uniX),length(uniColor),size(bTrace,2));

% fill response matrix
for iUni = 1:size(uni,1)
    indx = sum(allStims == repmat(uni(iUni,:),[size(allStims,1) 1]),2)==3;
    response(uni(iUni,1),uni(iUni,2),uni(iUni,3),:) = mean(bTrace(indx,:));
end

% get the trial averages
map = squeeze(response(:,:,2,:) + response(:,:,1,:));
map = permute(reshape(map,[size(map,1) size(map,2) sz(1) sz(2)]),[3 4 1 2]);
tuple.onpoff_rf = imresize(map,1/rz);
tuple.onpoff_rf = reshape(tuple.onpoff_rf,[size(tuple.onpoff_rf,1) size(tuple.onpoff_rf,2) size(map,3) size(map,4)]);
insert( obj, tuple );






