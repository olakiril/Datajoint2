function makeTuples( obj, key )

import vis2p.*

% get Optical data
disp 'loading movie...'
[name, date] = fetch1( Scans(key),'file_name','exp_date' );
if ~isempty(strfind(name,'.h5'));fend = [];else fend = '.h5';end
path = ['M:/IntrinsicImaging/' datestr(date, 'yymmdd')];
filename = getLocalPath([path '/' name fend]);
[data,Fs] = getOpticalData(filename);

disp 'trial separation...'
times = 1/Fs:1/Fs:size(data,1)/Fs;

%stim times and types
stim = fetch1(VisStims(key),'stim_file');

% Get stim data
startIndx = find(strcmp(vertcat(stim.eventTypes),'startTrial'));

% find trace segments
dataCell = cell(1,length(stim.params.trials));
for iTrial = 1:length(stim.params.trials)
    stimtime = stim.events(iTrial).syncedTimes(stim.events(iTrial).types == startIndx)/1000;
    dataCell{iTrial} = data(times>=stimtime & ...
        times<stimtime+ stim.params.trials(1).delayTime/1000,:,:);
end
data = dataCell;
clear dataCell

% remove incomplete trials
tracessize = cell2mat(cellfun(@size,data,'UniformOutput',0)');
indx = tracessize(:,1) >= stim.params.trials(1).delayTime/1000*9/10*Fs;
data = data(indx);
tracessize = tracessize(indx,1);

% equalize trial length
data = cellfun(@(x) permute(zscore(x(1:min(tracessize(:,1)),:,:)),[2 3 1]),data,'UniformOutput',0);
tf = Fs/size(data{1},3);
data = permute(cat(3,data{:}),[3 1 2]);

% get the vessel image
disp 'getting the vessels...'
k = [];
k.exp_date = key.exp_date;
k.mouse_id = key.mouse_id;
vesObj = Scans(k,'stim_engine = "None" and scan_prog = "Imager"');
if ~isempty(vesObj)
    scans = fetchn(vesObj,'scan_idx');
    [~,i] = min(abs(key.scan_idx - scans)); % select the closest scan
    keys = fetch(vesObj);
    name = fetch1(Scans(keys(i)),'file_name');
    vessels = squeeze(mean(getOpticalData([path '/' name '.h5'])));
else
    vessels = [];
end

imsize = size(data,2);
data = (bsxfun(@minus,data(:,:),mean(data(:,:)))); % subtract mean for the fft

T = 1/Fs; % frame length
L = size(data,1); % Length of signal
t = (0:L-1)*T; % time series

% get stimulation information
if strcmp(fetch1(VisStims(key),'exp_type'),'BarMappingExperiment')
    
%     stim = fetch1(VisStims(key),'stim_file');
%     s = ([stim.params.trials.sync]);
%     tf = 1/median(diff([s(:).start]))*1000;

    key.direction = stim.params.constants.trajectoryAngle;
else
    tf = 0.16;
    key.direction = 0;
end

% do it
disp 'computing...'
R = exp(2*pi*1i*t*tf)*data;
imP = squeeze(reshape((angle(R)),imsize,imsize));
imA = squeeze(reshape((abs(R)),imsize,imsize));

% save the data
disp 'inserting data...'
key.amp = imA;
key.ang = imP;
if ~isempty(vessels); key.vessels = vessels; end

insert(obj,key)

disp 'done!'

