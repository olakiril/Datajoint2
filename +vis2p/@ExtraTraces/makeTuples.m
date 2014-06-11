function makeTuples(obj,key)

import vis2p.*


% get tpReader
tp = tpReader(Scans(key));
if isempty(tp);return;end

% find eye file
eye_dir = findEye(obj,key);
if isempty(eye_dir) && strcmp(fetch1(Scans(key),'aim'),'shock') &&  getElectrodeChannelIndicies(tp)==2
    shock_trace = tp.elCh{2}(:);
    shock_trace = single(shock_trace>median(double(shock_trace))+std(double(shock_trace)));
    key.shock_trace = interp1(shock_trace,single(1:getSamplesPerFrame(getElectrodeChannel(tp,2)):length(shock_trace)));
    % insert data
    insert( obj, key );
% elseif ~isempty(eye_dir)
%     % Set parameters
%     pRad = [5 70]; % radii range for the pupil in pixels
%     inSens = 0.99;   % Sensitivity for circle detection
%     tic;
%     load([eye_dir '/eyeData.mat'])
%     eMask = eyeData.eMask;
%     xind = eyeData.xind:eyeData.xind+size(eMask,2)-1;
%     yind = eyeData.yind:eyeData.yind+size(eMask,1)-1;
%     save([eye_dir '/eyeData_old.mat'],'eyeData')
%     xyloObj = VideoReader([eye_dir '/eyemovie.avi']);
%     nFrames = xyloObj.NumberOfFrames;
%     eMask2 = normalize(imfilter(double(eMask),gausswin(250)*gausswin(250)','replicate'));
%     gfilt = gausswin(20)*gausswin(20)';
%     efilt = datenum(key.exp_date)<datenum('2013-09-20');
%     gfilt2 = gausswin(5)*gausswin(5)';
%     % initialize variables
%     Centers = nan(nFrames,2);
%     Radii = nan(nFrames,1);
%     
%     blobs = 1:ceil(nFrames/15):nFrames;
%     blobs(end+1) = nFrames+1;
%     for iblob = 1:length(blobs)-1;
%         frames = read(xyloObj,[blobs(iblob) blobs(iblob+1)-1]);
%         frames = single(squeeze(frames(yind,xind,1,:)));
%         centers = nan(size(frames,3),2);
%         radii = nan(size(frames,3),1);
%         isrt = blobs(iblob);
%         parfor indx = 1:size(frames,3);
%             if ~mod(indx,1000);display([num2str(indx+isrt) '/' num2str(nFrames)]);end
%             
%             % Read Frame
%             frame = frames(:,:,indx);
%             if efilt; frame = imfilter(frame+10,gfilt2);end
%             frame = 1./(frame.^2);
%             gframe = frame.*single(eMask2); % store this for the SNR computation
%             
%             % Filter & mask image
%             frame = normalize(imfilter(frame,gfilt,'replicate').*single(eMask2.^4));
%             
%             % Detect Edges
%             frame(frame<0.3) = 0.3;
%             frame(frame>0.7) = 0.7;
%             
%             % Detect Pupil with given range
%             [ce,r]=imfindcircles(frame,pRad,'method','TwoStage','sensitivity',inSens);
%             if isempty(r);[ce,r]=imfindcircles(frame,pRad,'sensitivity',1);end % increase threshold if nothing found
%             if isempty(r);centers(indx,:) = [nan nan];radii(indx)=nan;continue;end
%             
%             % Compute SnR of all the detected circles
%             SNR = nan(size(ce,1),1);
%             for ic = 1:size(ce,1)
%                 % pupil mask
%                 [x, y] = circlepoints(r(ic));
%                 pupil = roipoly(gframe,x+ce(ic,1),y+ce(ic,2));
%                 
%                 % surround
%                 [x, y] = circlepoints(r(ic)*2);
%                 sur = roipoly(gframe,x+ce(ic,1),y+ce(ic,2));sur(pupil)=false;
%                 
%                 % compute snr
%                 sig = mean(gframe(pupil));
%                 noise = mean(gframe(sur));
%                 SNR(ic) = sig/noise ;
%             end
%             
%             % sort by SNR
%             [~,isort] = sort(SNR,'descend');
%             ce = ce(isort,:);
%             r = r(isort);
%             
%             % store datapoints
%             centers(indx,:) = ce(1,:); % get the best SNR circle centers
%             radii(indx) =r(1); % get the best SNR circle radius
%             
%         end
%         Centers(blobs(iblob):blobs(iblob)+size(centers,1)-1,:) = centers;
%         Radii(blobs(iblob):blobs(iblob)+size(centers,1)-1) = radii;
%     end
%     % recreate the
%     disp 'Saving Data...'
%     eyeData.eMask = eMask;
%     eyeData.xind = xind(1);
%     eyeData.yind = yind(1);
%     eyeData.centers = Centers;
%     eyeData.radii = Radii;
%     eyeData.framerate = 1/(xyloObj.Duration/size(radii,1));
%     eyeData.duration = xyloObj.Duration;
%     save([eye_dir '/eyeData.mat'],'eyeData')
%     t = toc; fprintf('%s\n', ['. Elapsed time: ' num2str(round(t)) ' sec']);
    
elseif ~isempty(eye_dir)
    
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
    load([eye_dir 'eyeData.mat'])
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



