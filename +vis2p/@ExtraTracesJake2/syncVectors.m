function timeCorrector = syncVectors(obj,tpObj,H5File)

% function syncPhotodiodes(obj,trace1,trace2,fps1,fps2)
%
% Synchronizes 2 vector traces
%
% MF 2012-11-18

gains = -50:1:50; % range of gains in ms

% gain function
gainfix = @(x,gn,x1,xend)  (x-x1)*(1 + gn/(xend - x1)) + x1;

% load the tpfile related data
tpTS = double(readTimestamps(tpObj));
pd = getElectrodeChannel(tpObj,1);
pdSR = getSamplingRate(pd);
pdSpF = double(getSamplesPerFrame(pd));

% detect flip times
detflips = detectFlips(pd(:,1));
flips = interp1(tpTS, 1+(detflips/pdSpF),'linear','extrap'); 

% load the h5 file related data
tpTS =  loadHWS(H5File,'activity','times');
pd = loadHWS(H5File,'activity','photodiode');
pdSpF = length(pd)/length(tpTS);
pdSR = pdSpF/mean(diff(tpTS));
tpTS = tpTS*1000;
        
% detect swap times
detflips = detectFlips(pd);
swaps = interp1(tpTS, 1+(detflips/pdSpF),'linear','extrap');

% go through multiple gains 
amp = nan(length(gains),1); ind = amp;
for igain = 1:length(gains);
    
    % build stim vectors
    [vswaps, mdflips] = buildStimVector(gainfix(swaps,gains(igain),swaps(1),swaps(end)));
    vflips = buildStimVector(flips,mdflips);
    
    % remove trials that are not in both vectors from the end of the traces: 2012-09-07
    dvswaps = find(diff(vswaps)>0);
    dvflips = find(diff(vflips)>0);
    trind = min([length(dvswaps) length(dvflips)]);
    vswaps = vswaps(1:dvswaps(trind));
    vflips = vflips(1:dvflips(trind));
    
    % organize..
    trace{1} = vswaps;
    trace{2} = vflips;
    [maxL, maxI] = max([length(vswaps) length(vflips)]);
    [~, minI] = min([length(vswaps) length(vflips)]); 
    
    % find the optimal shift
    cor = xcorr(trace{maxI}+1,trace{minI}+1);
     
    [amp(igain), ind(igain)] = max(cor);

end
[~,icor] = max(amp);

delay = ind(icor) - maxL;
if maxI == 1;delay = -delay;end

% correct the times
timeCorrector = @(x) gainfix(x,gains(icor),swaps(1),swaps(end)) + ...
    flips(1) + delay -swaps(1);

% plot if nothing is requested
% if ~nargout
    figure
    plot(flips,ones(size(flips)),'.');
    hold on
    plot(timeCorrector(swaps),2*ones(size(swaps)),'.r')
    set(gca,'YLim',[0 20])
% end

% this function creates a time vector of 0s, and 1s where stimulus exists
function flips = detectFlips(flips)
flips = flips>mean(flips);
dflip = diff(flips,2);
flips(find(dflip==-2)+1) = 0;
flips(find(dflip==2)+1) = 1;
flips = find(flips);

function [stimVector, mdflips] = buildStimVector(stimTimes,mdflips)
cflips = round(stimTimes - stimTimes(1));
cflips(cflips == 0 ) = 10;
stimVector = ones(ceil(stimTimes(end) - stimTimes(1)),1)*10;
dflips = diff(cflips);
if nargin<2
    mdflips = median(dflips);
end
indx = find(dflips > 2*mdflips); % this assumes that off is 8x the on period
for istop = 1:length(indx)
    stimVector(cflips(indx(istop)):cflips(indx(istop)+1)) = 0;
end

