% Old Function
function timeCorrector = syncPhotodiodes(obj,tpObj,H5File,flipfreq)

% function syncPhotodiodes(obj,trace1,trace2,fps1,fps2)
%
% Synchronizes 2 Photodiode traces
%
% MF 2012-11-18
params.manual = 1;

if nargin<4
    flipfreq = 30; % flipping frequency of the stimulus.
end
gains = -1:.01:1; % range of gains in % of main clock

% gain function
gainfix = @(x,gn)  (x-x(1))*(1 + gn/100) + x(1);

% load the tpfile related data
tpTS = double(readTimestamps(tpObj));
pd = getElectrodeChannel(tpObj,1);
pdSR = getSamplingRate(pd);
pdSpF = double(getSamplesPerFrame(pd));

% detect flip times
detflips = detectFlipsM(pd(:,1),pdSR,flipfreq);
flips = interp1(tpTS, 1+(detflips/pdSpF),'linear','extrap');
ftime = tpTS(1);

% load the h5 file related data
tpTS =  loadHWS(H5File,'activity','times');
pd = loadHWS(H5File,'activity','photodiode');
pdSpF = length(pd)/length(tpTS);
pdSR = pdSpF/mean(diff(tpTS));
tpTS = tpTS*1000;

% detect swap times
detflips = detectFlipsM(pd,pdSR,flipfreq);
swaps = interp1(tpTS, 1+(detflips/pdSpF),'linear','extrap');

if params.manual
    [delay, gain]= alignVectors([flips ones(length(flips),1)],swaps+(ftime - swaps(1)),'marker','.');
    timeCorrector = @(x) gainfix(x,gain) + ftime + delay - swaps(1);
else
    % go through multiple gains
    amp = nan(length(gains),1); delay = amp;
    for igain = 1:length(gains);

        % build stim vectors
        [vswaps, mdflips] = buildStimVector(gainfix(swaps,gains(igain)));
        vflips = buildStimVector(flips,mdflips);

        % remove trials that are not in both vectors from the end of the traces: 2012-09-07
    %     dvswaps = find(diff(vswaps)>0);
    %     dvflips = find(diff(vflips)>0);
    %     trind = min([length(dvswaps) length(dvflips)]);
    %     vswaps = vswaps(1:dvswaps(trind));
    %     vflips = vflips(1:dvflips(trind));

        % organize..
        trace{1} = vswaps;
        trace{2} = vflips;
        [maxL, maxI] = max([length(vswaps) length(vflips)]);
        [~, minI] = min([length(vswaps) length(vflips)]);

        % find the optimal shift
        cor = xcorr(trace{maxI}+1,trace{minI}+1);
        [amp(igain), ind] = max(cor);
        delCor = [-1 1];
        delay(igain) = delCor(maxI)*(ind - maxL);
    end
    [~,icor] = max(amp);
    timeCorrector = @(x) gainfix(x,gains(icor)) + flips(1) + delay(icor) - swaps(1);
end
% plot if nothing is requested
% if ~nargout
    figure
    plot(flips,ones(size(flips)),'.');
    hold on
    plot(timeCorrector(swaps),1.5*ones(size(swaps)),'.r')
    set(gca,'YLim',[0 20])
    title(getFilename(tpObj),'interpreter','none')
% end

% this function creates a time vector of 0s, and 1s where stimulus exists
function [stimVector, mdflips] = buildStimVector(stimTimes,mdflips)
cflips = round(stimTimes - stimTimes(1));
cflips(cflips == 0 ) = 1;
stimVector = ones(ceil(stimTimes(end) - stimTimes(1)),1)*1;
dflips = diff(cflips);
if nargin<2
    mdflips = median(dflips);
end
indx = find(dflips > 8*mdflips); % this assumes that off is 8x the on period
for istop = 1:length(indx)
    stimVector(cflips(indx(istop)):cflips(indx(istop)+1)) = 0;
end

%% New function
% function timeCorrector = syncPhotodiodes(obj,tpObj,H5File,flipfreq,manual)
% 
% % function syncPhotodiodes(obj,trace1,trace2,fps1,fps2)
% %
% % Synchronizes 2 Photodiode traces
% %
% % MF 2012-11-18
% 
% if nargin<5
%     manual = 0;
% else
%     manual = 1;
%     marker = '.';
% end
% 
% if nargin<4
%     flipfreq = 30; % flipping frequency of the stimulus.
% end
% 
% % load the tpfile related data
% tpTS = double(readTimestamps(tpObj));
% pd = getElectrodeChannel(tpObj,1);
% pdSR = getSamplingRate(pd);
% pdSpF = double(getSamplesPerFrame(pd));
% 
% % detect flip times
% detflips = detectFlipsM(pd(:,1),pdSR,flipfreq);
% flips = interp1(tpTS, 1+(detflips/pdSpF),'linear','extrap');
% 
% % load the h5 file related data
% tpTS =  loadHWS(H5File,'activity','times');
% pd = loadHWS(H5File,'activity','photodiode');
% pdSpF = length(pd)/length(tpTS);
% pdSR = pdSpF/mean(diff(tpTS));
% tpTS = tpTS*1000;
% 
% % detect swap times
% detflips = detectFlipsM(pd,pdSR,flipfreq);
% swaps = interp1(tpTS, 1+(detflips/pdSpF),'linear','extrap');
% 
% if ~manual
%     % go through multiple gains
%     
%     % build stim vectors
%     [vswaps, mdflips] = buildStimVector(swaps);
%     vflips = buildStimVector(flips,mdflips);
%     
%     % remove trials that are not in both vectors from the end of the traces: 2012-09-07
%     %     dvswaps = find(diff(vswaps)>0);
%     %     dvflips = find(diff(vflips)>0);
%     %     trind = min([length(dvswaps) length(dvflips)]);
%     %     vswaps = vswaps(1:dvswaps(trind));
%     %     vflips = vflips(1:dvflips(trind));
%     
%     % organize..
%     trace{1} = vswaps;
%     trace{2} = vflips;
%     [maxL, maxI] = max([length(vswaps) length(vflips)]);
%     [minL, minI] = min([length(vswaps) length(vflips)]); %#ok<ASGLU>
%     
%     % find the optimal shift
%     cor = xcorr(trace{maxI},trace{minI});
%     [amp, ind] = max(cor); %#ok<ASGLU>
%     delay = ind - maxL;
%     % correct the times
%     timeCorrector = @(x)  x+ flips(1) + delay - swaps(1);
%     
% else
%     delay = alignVectors([flips ones(length(flips),1)],swaps +(tpTS(1) - swaps(1)),'marker',marker);
%     %         delay = alignData([times' rd],swaps +(tpTS(1) - swaps(1)),'marker',marker);
%     % correct the times
%     timeCorrector =  tpTS(1) + delay - swaps(1);
% end
% % plot if nothing is requested
% % if ~nargout
% figure
% plot(flips,ones(size(flips)),'.');
% hold on
% plot(timeCorrector(swaps),2*ones(size(swaps)),'.r')
% set(gca,'YLim',[0 20])
% % end
% 
% % this function creates a time vector of 0s, and 1s where stimulus exists
% function [stimVector, mdflips] = buildStimVector(stimTimes,mdflips)
% cflips = round(stimTimes - stimTimes(1));
% cflips(cflips == 0 ) = 10;
% stimVector = ones(ceil(stimTimes(end) - stimTimes(1)),1)*10;
% dflips = diff(cflips);
% if nargin<2
%     mdflips = median(dflips);
% end
% indx = find(dflips > 8*mdflips); % this assumes that off is 8x the on period
% for istop = 1:length(indx)
%     stimVector(cflips(indx(istop)):cflips(indx(istop)+1)) = 0;
% end
