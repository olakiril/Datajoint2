function makeTuples(obj,key)

% get tpReader
tp = tpReader(Scans(key));
if isempty(tp);return;end

% find eye file
eye_dir = findEye(obj,key);
if ~isempty(eye_dir)
    try
            load([eye_dir 'eyeDataJake.mat'])
    catch
        disp('Jake said no')
        return
    end
    
    fn = [eye_dir 'activity.h5'];
    
    try
        if ~isempty(VisStims(key))
            % fix timings from photodiode traces
            corfun = syncPhotodiodes(obj,tp,fn);
            timestamps = cell2mat( fetchn(VisStims(key),'frame_timestamps')')';
        else
            corfun = syncVectors(obj,tp,fn);
            timestamps = double(readTimestamps(tp));
        end
        
        % get some data
        times = loadHWS(fn,'activity','times')*1000;
        activity =loadHWS(fn,'activity','activity');
        activity(isinf(activity)) = nan;
        stimulation = loadHWS(fn,'activity','shock');
        
    catch
        disp('Corrupt activity file!')
        return
    end
    
    % mouse activity
    lact = length(activity);
    [times,indx]=unique(corfun(times));activity=activity(indx);
    key.speed_trace = interp1(times,activity,timestamps);
    
    % stimulations
    if ~isempty(stimulation)
        key.shock_trace = interp1(times,stimulation(indx),timestamps);
    else
        key.shock_trace = [];
    end
    
    % eye movements & pupil dilation
    load([eye_dir 'eyeDataJake.mat'])
    msort = @(x) permute(x,((size(x))==min(size(x)))+1);
    centers = msort(eyeData.centers);
    radii = msort(eyeData.radii);
    c1 = filterEye(centers(:,1));
    c2 = filterEye(centers(:,2));
    radii = filterEye(radii);
    xi = size(c1,1)/lact:size(c1,1)/lact:size(c1,1);
    centers =  interp1(1:size(c1,1),c1,xi,'spline')';
    centers(:,2)  =  interp1(1:size(c1,1),c2,xi,'spline');
    radii =  interp1(1:size(c1,1),radii,xi,'spline');
    key.eye_trace(:,1) = interp1(times,centers(indx,1),timestamps);
    key.eye_trace(:,2) = interp1(times,centers(indx,2),timestamps);
    key.pupil_trace = interp1(times,radii(indx),timestamps);
    key.quality = eyeData.quality;
        
    % whisker movements
    load([eye_dir 'whiskerData.mat'])
    whiskerData.shifts(1) = whiskerData.shifts(2); % first value is bad
    key.whisker_trace = interp1(times,whiskerData.shifts(indx),timestamps);
    
    % insert data
    insert( obj, key );
    
else
    return;% quit if no eye file is found
end

% filter eye movements
function ec = filterEye(x)

mx = nanmean(x);
x(isnan(x)) = mx;
x = x - mx;

% high amplitude noise
iters = 200;
thr =20; % max shift allowed between frames in pixels
for i = 1:iters
    d = [0; diff(x)];
    indb = ((d>thr&x>thr));
    inda = d<-thr&x<-thr ;
    indb(end) = false;indb(1) = false;
    inda(1) = false;inda(end) = false;
    x(inda) =x([inda(1);inda(1:end-1)]);
    x(indb) = x([indb(2:end);indb(end)]);
end

% low amplitude noise
thr = 3; % max shift allowed between frames in pixels
d = [0; diff(x)];
ind = d>thr&x>thr | d<-thr&x<-thr;
if sum(ind)>0;
    x(ind) = interp1(find(~ind),x(~ind),find(ind),'spline');
end

% filter
windowSize = 20;
xx = conv(x,ones(1,windowSize)/windowSize,'same');
xx(1:ceil(windowSize/2)) = mean(xx(1:ceil(windowSize/2)));
xx(end-ceil(windowSize/2):end) = mean(xx(end-ceil(windowSize/2):end));
ec = xx+mx;



