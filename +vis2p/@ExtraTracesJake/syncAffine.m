function times = syncAffine(obj,tpObj,H5File,flipfreq)

% function syncAffine(obj,trace1,trace2,fps1,fps2)
%
% Synchronizes 2 Photodiode traces
%
% MF 2013-11-18
params.manual = 0;

if nargin<4
    flipfreq = 30; % flipping frequency of the stimulus.
end

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
    % build stim vectors
    [vswaps, mdflips] = buildStimVector(swaps);
    vflips = buildStimVector(flips,mdflips);
    
    % organize..
    trace{1} = vswaps;
    trace{2} = vflips;
    [maxL, maxI] = max([length(vswaps) length(vflips)]);
    [~, minI] = min([length(vswaps) length(vflips)]);
    
    % find the optimal shift
    cor = xcorr(trace{maxI}+1,trace{minI}+1);
    [~, ind] = max(cor);
    delCor = [-1 1];
    delay = delCor(maxI)*(ind - maxL);
    timeCorrector = @(x) x + flips(1) + delay - swaps(1);
end

figure
plot(flips,ones(size(flips)),'.');
hold on
plot(timeCorrector(swaps),1.5*ones(size(swaps)),'.r')
set(gca,'YLim',[0 20])
title(getFilename(tpObj),'interpreter','none')

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

