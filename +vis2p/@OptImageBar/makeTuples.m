function makeTuples( obj, key )

import vis2p.*

vessels = [];

% get Optical data
disp 'loading movie...'
if strcmp(fetch1(Scans(key),'scan_prog'),'Imager')
    [name, date] = fetch1( Scans(key),'file_name','exp_date' );
    if ~isempty(strfind(name,'.h5'));fend = [];else fend = '.h5';end
    path = ['M:/IntrinsicImaging/' datestr(date, 'yymmdd')];
    filename = getLocalPath([path '/' name fend]);
    [data,Fs] = getOpticalData(filename);
    stim = fetch1(VisStims(key),'stim_file');
    times = 1/Fs:1/Fs:size(data,1)/Fs;
    
    % DF/F
    data(:,1,:) = data(:,2,:);
    %% process data



    mdata = mean(data);
    data = bsxfun(@rdivide,bsxfun(@minus,data,mdata),mdata);
%         data = permute(imresize(permute(data,[2 3 1]),1/4),[3 1 2]);
%     [c, p] = princomp(data(:,:));
%     data = p(:,2:end)*c(:,2:end)';
%     data = reshape(data,size(data,1),128,128);
else
    data = getFrames(Movies(key),1);
    vessels = mean(data,3);
    data = permute(data,[3 1 2]);
    Fs = fetch1(Movies(key),'fps');
    [stim, times] = fetch1(VisStims(key),'stim_file','frame_timestamps');
     times = times/1000;
     [path,name] = fetch1( Experiments*Scans(key), 'directory','file_name' );
    filename = getLocalPath([path '/' name]);
end


%%
%   traces =double(data(:,:));
%   
%         hp = 0.01;
% traces = traces + abs(min(traces(:)))+eps;
% traces = traces(:,:)./convmirr(traces(:,:),hamming(round(Fs/hp)*2+1)/sum(hamming(round(Fs/hp)*2+1)))-1;  %  dF/F where F is low pass
% traces = bsxfun(@plus,traces,abs(min(traces)))+eps;
% data = reshape(traces,size(data));

%%
disp 'trial separation...'

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

if strcmp(fetch1(Scans(key),'scan_prog'),'Imager')
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
    end
end

imsize = size(data);
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
elseif strcmp(fetch1(VisStims(key),'exp_type'),'FancyBarExperiment')
    dirs = [0 90];
    key.direction = dirs(stim.params.constants.Axis);
else
    tf = 0.16;
    key.direction = 0;
end

% do it
disp 'computing...'
R = exp(2*pi*1i*t*tf)*data;
imP = squeeze(reshape((angle(R)),imsize(2),imsize(3)));
imA = squeeze(reshape((abs(R)),imsize(2),imsize(3)));

% save the data
disp 'inserting data...'
key.amp = imA;
key.ang = imP;
if ~isempty(vessels); key.vessels = vessels; end

insert(obj,key)

disp 'done!'

