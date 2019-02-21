% use this function to detect the PSPs in your data using a PSP finding algorithm
% params = detectPSPs(data, PSPsDown, windows, windowSize, parameter pairs)
% or params = detectPSPs(data, PSPsDown, windows, windowSize)
% or params = detectPSPs(data, PSPsDown)
% if either pairs of parameter names and values xor a parameters structure is passed they must be the last arguement

% The algorithm used is based on that of Cohen and Miles 2000
% Briefly, the data are filtered and the derivative is taken.
% A function is generated based on the rule that a point is equal to the
% sum of the derivative at that point and the previous point only if the
% derivative is greater than (for upward PSPs) zero otherwise the function
% is zero at that point.
% This is then thresholded and sent for fitting by alpha functions and one
% variable is returned:
% Params(1,:) = amplitude of psp
% Params(2,:) = rise time
% Params(3,:) = x-offset of alpha function
% Params(4,:) = decay time

% To untangle multiple stacked PSPs as after a stim the readjustment of
% minYOffset/maxYOffset and also an iterative mechanism to break down giant

function [params decayFits] = detectPSPs(data, varargin)

warning off MATLAB:divideByZero
warning off all % to avoid a lot of errors with the best fitting

%%%
% Parse through the inputs
%%%

    % pick the orientation that we like
    if nargin > 0
        if size(data, 1) > size(data, 2)
            data = data';
        end
    else
    %     help detectPSPs
    %     error('Must specify data on which to operate');
        data = [];
    end

    % if more than one data set then reject
    if size(data, 1) > 1
        Disp('Input must be a 1 x N vector')
        return
    end

    % look for the structure of input parameters
    hasParameters = 0;
    for i = 1:nargin - 1
        if isstruct(varargin{i})
            hasParameters = i;
            parameters = varargin{i};
            break
        end
    end

    % find the first pair of input properties
    firstChar = nargin + 1;
    if ~hasParameters
        for i = 1:nargin - 1
            if ischar(varargin{i})
                firstChar = i + 1;
                break
            end
        end
    end

    switch firstChar - 1 - (hasParameters > 0)
        case 0
            windows = 1; % so that we include the whole time
            windowSize = length(data);
            % determine whether spikes are up or down
            if median(data) > mean(data)  
                %spikes are down
                PSPsDown = 1;
            else
                %spikes are up
                PSPsDown = 0;
            end
        case 1
            windows = 1; % so that we include the whole time
            windowSize = length(data);
            % determine whether spikes are up or down
            if median(data) > mean(data)
                %spikes are down
                PSPsDown = 1;
            else
                %spikes are up
                PSPsDown = 0;
            end
        case 2
            PSPsDown = varargin{1};
            windows = 1; % so that we include the whole time
            windowSize = length(data);
        case 3
            % determine whether spikes are up or down
            if median(data) > mean(data)
                %spikes are down
                PSPsDown = 1;
            else
                %spikes are up
                PSPsDown = 0;
            end
            windows = varargin{1};
            windowSize = varargin{2};        
        case 4
            if ismember(varargin{1}, [0 1])
                PSPsDown = varargin{1};
                windows = varargin{2};
                windowSize = varargin{3};
            else
                % determine whether spikes are up or down
                if median(data) > mean(data)
                    %spikes are down
                    PSPsDown = 1;
                else
                    %spikes are up
                    PSPsDown = 0;
                end
                windows = varargin{1};
                windowSize = varargin{2};  
                parameters = varargin{3};            
            end
        case 5
            PSPsDown = varargin{1};
            windows = varargin{2};
            windowSize = varargin{3};
            parameters = varargin{4};
    end

    if ~hasParameters 
        if PSPsDown
            parameters = struct('minAmp', (-5000),... % minimum allowable amplitude for alpha functions (in units of samples)
                'maxAmp', (-5),... % maximum allowable amplitude
                'minTau', (5),... % minimum allowable tau for alpha functions (in units of samples)
                'maxTau', (100),... %maximum allowable tau
                'minYOffset', (-100),... % minimum allowable yOffset for alpha functions (in units of mV)
                'maxYOffset', (-30),... % maximum allowable yOffset
                'minDecay', (5),... % minimum allowable decay tau
                'maxDecay', (500),... % maximum allowable decay tau
                'derThresh', (1),... % threshold used to determine if the change of derivative is great enough to separate out two alpha functions
                'closestEPSPs', (5),... % second EPSP is not recognized if it is closer than this value to the previous one (in units of samples)
                'errThresh', (0.08),... % threshold for standard error above which a fit with multiple alphas is attempted
                'dataFilterType', 1,... % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                'derFilterType', 3,... % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                'dataFilterLength', 5,... % length of data filter
                'derFilterLength', 7,... % length of derivative filter
                'debugging', (0),... % if set to 1 then debugging figures appear
                'dataStart', 1,... % index of first data point
                'forceDisplay', 0,... % forces a graphical output even if other outputs are taken
                'noFit', 0); % turns off best fitting of start time and decay Tau when 1
        else
            parameters = struct('minAmp', (5),... % minimum allowable amplitude for alpha functions (in units of samples)
                'maxAmp', (5000),... % maximum allowable amplitude
                'minTau', (5),... % minimum allowable tau for alpha functions (in units of samples)
                'maxTau', (100),... %maximum allowable tau
                'minYOffset', (-100),... % minimum allowable yOffset for alpha functions (in units of mV)
                'maxYOffset', (-30),... % maximum allowable yOffset
                'minDecay', (5),... % minimum allowable decay tau
                'maxDecay', (500),... % maximum allowable decay tau            
                'derThresh', (1),... % threshold used to determine if the change of derivative is great enough to separate out two alpha functions
                'closestEPSPs', (5),... % second EPSP is not recognized if it is closer than this value to the previous one (in units of samples)
                'errThresh', (0.004),... % threshold for standard error above which a fit with multiple alphas is attempted
                'dataFilterType', 1,... % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                'derFilterType', 3,... % 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                'dataFilterLength', 5,... % length of data filter
                'derFilterLength', 7,... % length of derivative filter
                'debugging', (0),... % if set to 1 then debugging figures appear
                'dataStart', 1,... % index of first data point
                'forceDisplay', 0,... % forces a graphical output even if other outputs are taken
                'noFit', 0); % turns off best fitting of start time and decay Tau when 1
        end
    end

    if firstChar < nargin
        if ~mod(nargin - firstChar, 2)
            disp('Inappropriate input');
            return
        else
            % specific parameters where passed so figure out which they were
            strArray = fieldnames(parameters);
            for paramIndex = firstChar - 1:2:numel(varargin)
                whichParam = strfind(strArray, varargin{paramIndex});
                if find(cell2mat(whichParam))
                    parameters.(varargin{paramIndex}) = varargin{paramIndex + 1};
                else
                    error([varargin{paramIndex} ' is not a valid parameter'])
                end
            end
        end
    end

% make sure that the filter length is odd
if parameters.dataFilterLength / 2 == fix(parameters.dataFilterLength / 2)
    parameters.dataFilterLength = parameters.dataFilterLength + 1;
end

if parameters.derFilterLength / 2 == fix(parameters.derFilterLength / 2)
    parameters.derFilterLength = parameters.derFilterLength + 1;
end

% reset the error threshold to reflect the noise in the signal
sortedData = sort(data);
if PSPsDown
    dataMean = mean(sortedData(size(data, 2) * .7:size(data, 2) * .9)); % mean of top 10-30% of data
    dataDev = std(sortedData(size(data, 2) * .6:size(data, 2) * .8)); % standard deviation of top 20-40% of data
    % or noise = data - dataFilt; dataDev = std(noise(1:1000));
    % dataDev = std(highPass(data)) * 3;
else
    dataMean = mean(sortedData(size(data, 2) * .1:size(data, 2) * .3)); % mean of bottom 10-30% of data
    dataDev = std(sortedData(size(data, 2) * .2:size(data, 2) * .4)); % standard deviation of bottom 20-40% of data
    %dataDev = std(highPass(data)) * 3;
end
parameters.errThresh = parameters.errThresh / dataDev;

% determine if data is current or voltage clamp
if nargin < 5 && ~hasParameters
    if max(data) > 50 || min(data) < -100
        % this must be voltage clamp
        parameters.minYOffset = dataMean - 200 * dataDev;
        parameters.maxYOffset = dataMean + 200 * dataDev;
    else
        % probably current clamp
        parameters.derThresh = .3;
        parameters.minAmp = 0.2;
        parameters.maxAmp = 50;
        parameters.minYOffset = dataMean - 3;
        parameters.maxYOffset = dataMean + 3;
    end
end

% if the person mixed up the max and min values then switch them
if (PSPsDown && (parameters.minAmp > 0 || parameters.maxAmp > 0)) || (~PSPsDown && (parameters.minAmp < 0 || parameters.maxAmp < 0))
    tempAmp = min(sign(-PSPsDown + 0.5) * abs([parameters.minAmp parameters.maxAmp]));
    parameters.maxAmp = max(sign(-PSPsDown + 0.5) * abs([parameters.minAmp parameters.maxAmp]));
    parameters.minAmp = tempAmp;
end

% filter the raw data
switch parameters.dataFilterType
    case 0 % no filtering
        dataFilt = data;
    case 1 % window filtering
        dataFilt = movingAverage(data, parameters.dataFilterLength);
    case 2 % median filtering
        dataFilt = medfilt1(data, parameters.dataFilterLength);
    case 3 % savitsky-golay
        dataFilt = sgolayfilt(data, 2, parameters.dataFilterLength);
end

% filter the derivative
tempDiff = diff(dataFilt);
switch parameters.derFilterType
    case 0 % no filtering
        dataDerFilt = tempDiff;
    case 1 % window filtering
        dataDerFilt = movingAverage(tempDiff, parameters.derFilterLength);
    case 2 % median filtering
        dataDerFilt = medfilt1(tempDiff, parameters.derFilterLength);
    case 3 % savitsky-golay
        dataDerFilt = sgolayfilt(tempDiff, 2, parameters.derFilterLength);
end
clear tempDiff;

% filter as per Cohen and Miles 2000
outData = zeros(size(dataFilt));
if PSPsDown
    for index = 2:length(dataFilt)
        if dataDerFilt(index - 1) < 0
            outData(index) = outData(index - 1) + dataDerFilt(index - 1);
        end
    end
else
    for index = 2:length(dataFilt)
        if dataDerFilt(index - 1) > 0
            outData(index) = outData(index - 1) + dataDerFilt(index - 1);
        end
    end    
end

if PSPsDown
    % find where derivative of this function is changing from negative to positive
    functionDer = diff(outData);
    peaks = find((functionDer(2:length(functionDer)) ./ functionDer(1:length(functionDer) -1) < 0 | functionDer(2:length(functionDer)) == 0) & functionDer(1:length(functionDer) - 1) < 0);
else
    % find where derivative of this function is changing from positive to negative
    functionDer = diff(outData);
    peaks = find((functionDer(2:length(functionDer)) ./ functionDer(1:length(functionDer) -1) < 0 | functionDer(2:length(functionDer)) == 0) & functionDer(1:length(functionDer) - 1) > 0);
end

% for each such value greater than derThresh find where the function last
% began to deviate from 0 and call that an event start
numStarts = 0;
whereNull = find(outData == 0);
whereStarts = zeros(length(peaks), 1);
wherePeaks = whereStarts;
for index = 1:length(peaks)
    if abs(outData(peaks(index))) > parameters.derThresh
        numStarts = numStarts + 1;
        whereStarts(numStarts) = whereNull(find(whereNull < peaks(index), 1, 'last'));
        wherePeaks(numStarts) = peaks(index);
    end
end
whereStarts(numStarts + 1) = length(outData);
whereStarts(numStarts + 2:end) = [];
wherePeaks(numStarts + 1:end) = [];

if parameters.debugging
    figure
    plot(data*0.4)
    hold on
    plot(dataFilt);
    plot(dataDerFilt);
    line(1:length(outData), outData + mean(dataFilt), 'Color', [1 0 0]);
    line(1:length(functionDer), functionDer + mean(dataFilt), 'Color', [0 1 0]);
    line(wherePeaks, dataFilt(wherePeaks), 'Color', [1 0 1], 'linestyle', 'none', 'marker', '+');
    line(whereStarts, dataFilt(whereStarts), 'Color', [0 0 0], 'linestyle', 'none', 'marker', '+');
    if PSPsDown
        line([1 length(dataFilt)], -[parameters.derThresh parameters.derThresh] + mean(dataFilt), 'Color', [1 0 0], 'linestyle', ':');
    else
        line([1 length(dataFilt)], [parameters.derThresh parameters.derThresh] + mean(dataFilt), 'Color', [1 0 0], 'linestyle', ':');        
    end
    fitLine = line([0 1], [dataMean dataMean], 'Color', [0 0 0]);
    decayLine = line([0 1], [dataMean dataMean], 'Color', [1 0 1]);
    riseStart = line([0 1], [dataMean dataMean], 'linestyle', 'none', 'marker', 'x', 'Color', [0.5430 0.2695 0.0742]);
    fitProps = annotation('textbox',[.8 .1 .1 .1], 'backgroundcolor', [1 1 1]);
    legend('Data', 'Up-Only Function', 'Derivative of U-O Fn', 'Event Peaks', 'Events Starts', 'Threshold for choosing', 'Fit line', 'Decay Fit', 'Rise Line');
end

% set a counter for events
numEvents = 0;
% create output array for alpha best fit parameters
params = zeros(1,4);

for indexEPSP = 1:length(whereStarts) - 1
    if whereStarts(indexEPSP + 1) - whereStarts(indexEPSP) < parameters.closestEPSPs
        continue % skip this iteration if the next fitting spot is really close
    end
    
    if any(~(whereStarts(indexEPSP)  < windows  - windowSize | whereStarts(indexEPSP)  > windows +  2 * windowSize)) % skip this iteration if the found event is outside of the given windows + a buffer of windowSize on each side in case the location is imprecise
        % must get tau within +/- 50% of the true value or the fitting
        % algorithm does not converge (norm of step goes below boundary
        % condition)
        tauGuess = max(wherePeaks(indexEPSP) - whereStarts(indexEPSP), parameters.minTau);
        ampGuess = dataFilt(wherePeaks(indexEPSP)) - dataFilt(whereStarts(indexEPSP));
        whereStop = whereStarts(indexEPSP + 1) - whereStarts(indexEPSP);
        if whereStop > 4 * tauGuess
            whereStop = 2 * tauGuess;
        end

        [bestParams, residuals] = nlinfit(whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop, dataFilt(whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop), @alpha, [ampGuess, tauGuess, whereStarts(indexEPSP), dataFilt(whereStarts(indexEPSP))], statset('MaxIter', 1000, 'FunValCheck', 'off'));
        %calculate standard error
        bestStdErr = sqrt(residuals' * residuals) / (length(residuals) - 2) / abs(bestParams(1));

        %show the fit line
        if parameters.debugging
            set(fitLine, 'Xdata', whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop, 'Ydata', alpha(bestParams, whereStarts(indexEPSP):whereStarts(indexEPSP) + whereStop));
            xBounds = get(gca, 'Xlim');
            if xBounds(2) < whereStarts(indexEPSP) + whereStop
               set(gca, 'Xlim', [whereStarts(indexEPSP) - 50 whereStarts(indexEPSP) + 250]);
            end
        end        
        
        if ~(parameters.minAmp < bestParams(1)  && parameters.maxAmp * 2.7183 > bestParams(1)  &&... % is the amplitude ok?
                    parameters.minTau < bestParams(2)  && parameters.maxTau > bestParams(2)  &&... % is the tau ok?
                    whereStarts(indexEPSP + 1) - bestParams(3) > bestParams(2) * 0.5 &&... % do we really have enough to fit a function to it?
                    parameters.minYOffset < bestParams(4) && parameters.maxYOffset > bestParams(4)) ||...% is it in the right place? 
                    ~any(~(bestParams(3)  < windows  | bestParams(3)  > windows + windowSize)) % is it in the window that we're looking in
            bestStdErr = inf; % set high so that any other fit is accepted since this one has parameters out of range
        end
        
        % refit variables
        if bestStdErr < parameters.errThresh       
            
            bestParams(1) = dataFilt(wherePeaks(indexEPSP)) - dataFilt(whereStarts(indexEPSP));  
            bestParams(4) = dataFilt(whereStarts(indexEPSP));
            bestParams(2) = wherePeaks(indexEPSP) - bestParams(3); %whereStarts(indexEPSP);
            if ~parameters.noFit 
                % change bestParams(4) to the decay tau
                whichData = round(bestParams(3) + bestParams(2) * 1.5:min([bestParams(3) + bestParams(2) * 1.5 + 100 whereStarts(indexEPSP + 1)]));
                if length(whichData) > 5
                    [bestParams(4) FittedDecay] = fitDecaySingle(dataFilt(whichData));
                    residuals = FittedDecay' - dataFilt(whichData);
                    if sqrt(residuals' * residuals) / (length(whichData) - 2) / abs(bestParams(1)) > parameters.errThresh
                        bestParams(4) = NaN;
                        FittedDecay = NaN;
                    end
                else
                    bestParams(4) = NaN;
                    FittedDecay = NaN;
                end

                if isnan(bestParams(4))
                   % make a manual guess
                   if PSPsDown == 0
                        tempDecay = find(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)])) < (dataFilt(wherePeaks(indexEPSP)) - log(2) * (dataFilt(wherePeaks(indexEPSP)) - min(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)]))))), 1, 'first');
                   else
                        tempDecay = find(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)])) > (dataFilt(wherePeaks(indexEPSP)) - log(2) * (dataFilt(wherePeaks(indexEPSP)) - max(dataFilt(wherePeaks(indexEPSP) + 2:min([wherePeaks(indexEPSP) + 90 whereStarts(indexEPSP + 1)]))))), 1, 'first');                       
                   end
                   if ~isempty(tempDecay)
                        bestParams(4) = tempDecay;
                        FittedDecay = ones([1 length(whichData)]) * dataFilt(wherePeaks(indexEPSP));
                   else
                        bestParams(4) = parameters.minDecay - 1;
                   end
                end

                if parameters.debugging               
                   set(fitProps, 'String',  [{['\bf Amp: ' num2str(bestParams(1))]};{[' Rise: ' num2str(bestParams(2))]};{[' Start: ' num2str(bestParams(3))]};{[' Decay: ' num2str(bestParams(4))]}]);
                   set(decayLine, 'Xdata', whichData, 'Ydata', FittedDecay);
                   set(gca, 'ylimmode', 'auto');
                   set(gca, 'xlim', [max([1 whereStarts(indexEPSP) - 10]) min([max([1 whereStarts(indexEPSP) - 10]) + max([bestParams(2) * 5 200]) length(data)])]);
                   savedLims = get(gca, 'ylim');
                   set(riseStart, 'xdata', bestParams(3), 'ydata', dataFilt(round(bestParams(3))));
                   set(gca, 'ylim', savedLims);
                   pause
                end
            else
                bestParams(4) = parameters.minDecay + 1;
            end
            
            if bestParams(4) >= parameters.minDecay && bestParams(4) <= parameters.maxDecay &&...
                    bestParams(1) >= parameters.minAmp && bestParams(1) <= parameters.maxAmp
                params(numEvents + 1, :) = bestParams;
                if ~parameters.noFit
                    decayFits{numEvents + 1} = FittedDecay;
                end
                numEvents = numEvents + 1;
            end
        end
    end
end

if parameters.noFit
    decayFits = nan(size(params));
    params(:,4) = nan;
end

% show how well our found trace fits the original data
if nargout == 0 || parameters.forceDisplay
    if max(data) > 50 || min(data) < -100
        parentHandle = scope(data, 1:length(data));
    else
        parentHandle = scope(data, 1:length(data));
    end

    set(0, 'currentfigure', get(parentHandle, 'parent'));
    set(gcf, 'currentAxes', parentHandle);      
    
    try
        for i = 1:size(params, 1)
            % plot starts
            line(parameters.dataStart + params(i, 3), data(round(params(i, 3))), 'Color', [0 0 1], 'linestyle', 'none', 'marker', '+', 'buttondownfcn', ['set(ancestor(gcbo, ''figure''), ''name'', ''Start time = ' sprintf('%6.1f', params(i,3)) ''')']);
            % plot amplitudes
            line(parameters.dataStart + [params(i, 3) + params(i, 2) params(i, 3) + params(i, 2)], [data(round(params(i, 3) + params(i, 2))) - params(i, 1) data(round(params(i, 3) + params(i, 2)))], 'Color', [0 1 0], 'buttondownfcn', ['set(ancestor(gcbo, ''figure''), ''name'', [''Rise time = ' sprintf('%4.1f', params(i,2)) ', Amplitude = ' sprintf('%7.2f', params(i,1)) '''])']);
            if ~parameters.noFit
                % plot decay tau
                line(parameters.dataStart + ((params(i, 3) + params(i, 2) * 1.5:params(i, 3) + params(i, 2) * 1.5 + numel(decayFits{i}) - 1)), decayFits{i}, 'Color', [1 0 0], 'buttondownfcn', ['set(ancestor(gcbo, ''figure''), ''name'', ''Decay tau = ' sprintf('%6.2f', params(i,4)) ''')']);        
            end
        end
    catch
        % at least return something
        
    end
end
end

function outData = movingAverage(inData, windowSize)
% use a boxcar filter of length windowSize points on inData
% filteredData = movingAverage(rawData, windowSize);
% defaults:
%   windowSize = 10 points

if nargin < 2
    windowSize = 10;
end

if size(inData, 1) > size(inData, 2)
    longSide = 1;
    flatData = ones(windowSize, 1);
else
    longSide = 2;
    flatData = ones(1, windowSize);
end

cheatShift = int32(windowSize / 2);
outData = filter(flatData./(windowSize),1,cat(longSide, flatData.*inData(1), inData, flatData.*inData(end)));
outData = outData(windowSize + cheatShift:length(inData) + windowSize + cheatShift - 1);
end

% function output = alpha([Amp, tau, xOffset, yOffset], xData)
%
%    (t-xOffset)     -(t-xOffset)
% A* ----------- * e^ ---------- + yOffset
%       tau               tau

function output = alpha(Params, x)

% see if we were passed an array of alphas or just one
if size(Params, 1) > 1
    output = zeros(1, length(x)) + mean(Params(:,4));
    yOffset = 0;
    for index = 1:size(Params, 1)
        Amp = Params(index,1);
        tau = Params(index,2);
        xOffset = Params(index,3);
        
        %eliminate pesky errors
        if tau == 0
            tau = 0.00000000001;
        end        
        tempOutput = Amp * (x - xOffset) / tau .* exp(-(x - xOffset) / tau) + yOffset;

        % make any x parts of the alpha function on the opposite side of the baseline from the alpha equal to the baseline
        if Amp > 0
            tempOutput(tempOutput < yOffset) = yOffset;
        else
            tempOutput(tempOutput > yOffset) = yOffset;
        end
        output = output + tempOutput;
    end
else
    Amp = Params(1);
    tau = Params(2);
    xOffset = Params(3);
    yOffset = Params(4);

        %eliminate pesky errors
        if tau == 0
            tau = 0.00000000001;
        end
        output = Amp * (x - xOffset) / tau .* exp(-(x - xOffset) / tau) + yOffset;

        % make any x parts of the alpha function on the opposite side of the baseline from the alpha equal to the baseline
        if Amp > 0
            output(output < yOffset) = yOffset;
        else
            output(output > yOffset) = yOffset;
        end
end
end

function [decayTau1 FittedCurve estimates sse] = fitDecaySingle(yData, PSPtype)
% fits tau to ydata that is slowing in x (only determined by the initial offset)
% [decayTau1 FittedCurve estimates] = fitDecaySingle(yData, PSPtype);

if nargin == 1
   if yData(end) >  yData(1)
       PSPtype = -1;
   else
       PSPtype = 1;
   end
end

if size(yData, 1) < size(yData, 2)
    yData = yData';
end

originalLength = length(yData);
if length(yData) > 100000
    yData(fix(length(yData)/ 500) * 500 + 1:end) = [];
    yData = mean(reshape(yData, 500, []))';
    downSampling = 500;
elseif length(yData) > 10000
    yData(fix(length(yData)/ 50) * 50 + 1:end) = [];                
    yData = mean(reshape(yData, 50, []))';
    downSampling = 50;
else
    downSampling = 1;
end

if length(yData) < 6
    decayTau1 = NaN;
    FittedCurve = NaN;
    estimates = NaN;
    sse = NaN;
    return
end

% generate xData
xdata = (0:length(yData) - 1)';

% Call fminsearch with a random starting point.
start_point = [yData(end) yData(1) - yData(end) length(yData) * -.7];
model = @expfun;
estimates = fminsearch(model, start_point, optimset('MaxFunEvals', 1000, 'Display', 'none'));

% check these for realism

if length(yData) < 200 && (abs(estimates(1) - yData(end)) > abs(yData(1) / 2) || (((estimates(3) >= 0) || estimates(2) * PSPtype == -1)))
   % try fitting with a single exponential
   % [decayTau1 FittedCurve] = fitDecayShort(yData, PSPtype);
   
   % try shaving a few points off of the ends
   if length(yData) > 7
       [decayTau1 FittedCurve estimates sse] = fitDecaySingle(yData(2:end - 2), PSPtype);
       %estimates(1) = estimates(1) + xdata(2) - xdata(5);
       if isnan(decayTau1)
           return
       end
   else
       decayTau1 = NaN;
       FittedCurve = NaN;
       estimates = NaN;
       sse = NaN;
       return
   end
else
    % return exponent with greatest magnitude
    decayTau1 = abs(estimates(3));
end

decayTau1 = decayTau1 * downSampling;
estimates(3) = estimates(3) * downSampling;

if nargout == 0
    figure, plot(xdata .* PSPtype, yData)
    if isnan(decayTau1)
        annotation('textbox',[.25 .5 .5 .1], 'backgroundcolor', [1 1 1], 'String', 'Unable to fit', 'edgecolor', 'none', 'fontsize', 24, 'horizontalalignment', 'center', 'verticalalignment', 'middle');
    else
        line(xdata .* PSPtype, estimates(1) + PSPtype * estimates(2) .* exp(xdata./estimates(3)), 'color', [1 0 0]);
        set(gcf, 'numbertitle', 'off', 'name', ['Tau = ' sprintf('%4.2f', decayTau1)]);    
    end
elseif nargout > 1
    xdata = (0:originalLength - 1)';
    FittedCurve = estimates(1) + PSPtype * estimates(2) .* exp(xdata./estimates(3));
end

if nargout > 3
    if downSampling == 1
        sse = expfun(estimates);
    else
        sse = nan;
    end
end

% expfun accepts curve parameters as inputs, and outputs sse,
% the sum of squares error for: offest + A * exp(-xdata./tau1) + B * exp(-xdata./tau2) - yData
    function sse = expfun(params)
        ErrorVector = (params(1) + PSPtype * params(2) .* exp(xdata./params(3))) - yData;
        sse = ErrorVector' * ErrorVector;
    end
end

function handles = scope(name, data, xData)
% display data channels in a scrollable window
% scope(chanData)
% scope(name, chanData)
% scope(chanData, timeData)
% scope(name, chanData, timeData)

% left click and drag on the axes to zoom in, right click to zoom out.
% clicking on the display will set or unset the cursors.
% the text box of the bottom left corner is the zoom factor

% parse input
switch nargin
    case 1
        data = name;
        name = 'Scope';
        xData = 1:length(data);
    case 2
        if ischar(name)
            xData = 1:length(data);
        else
            xData = data;
            data = name;
            name = 'Scope';
        end
    case 3
        % do nothing
    otherwise
        disp('Improper input')
end

if size(data, 1) < size(data, 2)
    data = data';
end

if size(data) == 1
    noData = 1;
    data = ones(1,data);
else
    noData = 0;
end

%initialize plotting window
set(0,'Units','normal');

figure('NumberTitle','off',...
    'Name', name,...
    'menu', 'none',...
    'Units','normal',...
    'Position',[0 .025 1 .95],...
    'Visible', 'on',...
    'UserData', -100,...
    'windowButtonDownFcn', @mouseDownScope,...
    'WindowButtonMotionFcn', @movePointers,...
    'windowButtonUpFcn', @mouseUpScope);
clear control_info;

uicontrol('Style','slider','Units','normal','Position', [0 0 1 .015], 'value', 0, 'sliderStep', [1 1/0], 'callback', @traceScroll);
uicontrol('Style','edit','Units','normal','Position',[0 .018 .02 .02], 'string', '1', 'callback', @traceScroll);

for index = 1:size(data, 2);
    handles(index) = axes('Position', [.05 .05 + (index - 1) / (size(data, 2) / .95)  .95 .95 / size(data, 2)], 'nextplot', 'add');
    startLine = line([1 1], [0 1], 'color', [0 0 0], 'userData', 0);
    stopLine = line([0 0], [0 1], 'color', [0 0 0]);
    lineText = text(1, 1, '0', 'color', [0 0 0], 'VerticalAlignment', 'bottom', 'linestyle', 'none');
    text(0, 1, '0', 'visible', 'off', 'VerticalAlignment', 'bottom', 'color', [0 0 0], 'linestyle', 'none')
    if index > 1
        set(gca, 'xticklabel', '', 'xtick', []);
    end
    if ~noData
        line(xData, data(:, size(data, 2) + 1 - index));
        if max(data(:, size(data, 2) + 1 - index)) > min(data(:, size(data, 2) + 1 - index))
            set(gca, 'ylim', [min(data(:, size(data, 2) + 1 - index)) max(data(:, size(data, 2) + 1 - index))]);
        else
            set(gca, 'ylim', [data(1, size(data, 2) + 1 - index) - 1 data(1, size(data, 2) + 1 - index) + 1]);
        end
        set(gca, 'xlim', [min(xData) max(xData)], 'userData', lineText);
        set(startLine, 'yData', get(gca, 'ylim'));
        set(stopLine, 'yData', get(gca, 'ylim'));
    end
end

set(gcf, 'userData', 0);
end

function mouseDownScope(varargin)
    pointerLoc = get(gcf, 'CurrentPoint');
    imageLoc = get(gca, 'Position');
    if pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) > .05
        switch get(gcf, 'SelectionType')
            case 'normal' %left mouse button clicked
                kids = get(gcf, 'children');
                firstKid = get(kids(1), 'children');
                if get(firstKid(length(firstKid)), 'userData') == 0
                    for index = 1:length(kids) - 2
                        kidLines = get(kids(index), 'children');
                        set(kidLines(length(kidLines)), 'color', [1 0 0]);
                    end
                    set(firstKid(length(firstKid)), 'userData', 1);   
                else
                    for index = 1:length(kids) - 2
                        kidLines = get(kids(index), 'children');
                        set(kidLines(length(kidLines)), 'color', [0 0 0]);
                        set(kidLines(length(kidLines) - 1), 'xData', [min(get(kids(index), 'xlim')) min(get(kids(index), 'xlim'))]);
                        set(kidLines(length(kidLines) - 3), 'visible', 'off')
                    end
                    set(firstKid(length(firstKid)), 'userData', 0);                
                end
            case 'extend' %middle mouse button clicked

            case 'alt' %right mouse button clicked

            case 'open' %double click

        end
    elseif pointerLoc(1) <= imageLoc(1) && pointerLoc(2) > 0.05
        kids = get(gcf, 'children');
        whichAxis = fix((pointerLoc(2) - .05) / .95 * (length(kids) - 2)) + 1;
        switch get(gcf, 'SelectionType')
            case 'normal' %left mouse button clicked
                set(gcf, 'userData', [whichAxis pointerLoc(2)]);
                uicontrol('units', 'normal', 'Style', 'text', 'String', '', 'backgroundColor', [0 0 1], 'Position', [0, pointerLoc(2), .01, 0.0001]);
            case 'extend' %middle mouse button clicked

            case 'alt' %right mouse button clicked

            case 'open' %double click

        end
    elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) <= .05
        %kids = get(gcf, 'children');
        switch get(gcf, 'SelectionType')
            case 'normal' %left mouse button clicked
                set(gcf, 'userData', [-1 pointerLoc(1)]);
                uicontrol('units', 'normal', 'Style', 'text', 'String', '', 'backgroundColor', [0 0 1], 'Position', [pointerLoc(1), .015, 0.0001, .01]);
            case 'extend' %middle mouse button clicked

            case 'alt' %right mouse button clicked

            case 'open' %double click

        end
    end
end

function movePointers(varargin)
    if get(gcf, 'userData') ~= -100
        pointerLoc = get(gcf, 'CurrentPoint');
        imageLoc = get(gca, 'Position');
        if length(get(gcf, 'userData')) > 1    
            info = get(gcf, 'userData');
            blueRect = get(gcf, 'children');
            blueRect = blueRect(1);
            if info(1) == -1
                if pointerLoc(1) > info(2)
                    set(blueRect, 'position', [info(2), .015, pointerLoc(1) - info(2), .01]);
                elseif pointerLoc(1) < info(2)
                    set(blueRect, 'position', [pointerLoc(1), .015, info(2) - pointerLoc(1), .01]);
                end
            else
                if pointerLoc(2) > info(2)
                    set(blueRect, 'position', [0, info(2), .01, pointerLoc(2) - info(2)]);
                elseif pointerLoc(2) < info(2)
                    set(blueRect, 'position', [0, pointerLoc(2), .01, info(2) - pointerLoc(2)]);
                end
            end
        elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) > .05
            set(gcf, 'pointer', 'crosshair');
            kids = get(gcf, 'children');
            firstKid = get(kids(1), 'children');
            whichAxis = fix((pointerLoc(2) - .05) / .95 * (length(kids) - 2)) + 1;
            axisKids = get(kids(length(kids) - 1 - whichAxis), 'children');
            xCoord = (pointerLoc(1) - .05) / .95  * diff(get(gca,'Xlim')) + min(get(gca,'Xlim'));
            xData = get(axisKids(end - 4), 'xData');
            [junk whereX] = min(abs(xCoord - xData));
            
           % yCoord = min(abs(yCoords(xCoord) - round(((pointerLoc(2) - .05 - (whichAxis - 1) * (.95 / (length(kids) - 2))) / (.95 / (length(kids) - 2)))  * diff(get(kids(length(kids) - 1 - whichAxis),'ylim')) + min(get(kids(length(kids) - 1 - whichAxis),'ylim')))));     
            
            if get(firstKid(length(firstKid)), 'userData')   == 0
                for index = 1:length(kids) - 2
                    kidLines = get(kids(index), 'children');
                    set(kidLines(length(kidLines)), 'xData', [xData(whereX) xData(whereX)], 'yData', get(kids(index), 'ylim'));
                    yData = get(kidLines(length(kidLines) - 4), 'yData');
                    yBounds = get(kids(index), 'ylim');
                    if abs(yBounds(1) - yData(whereX)) > abs(yBounds(2) - yData(whereX))
                        whereY = yBounds(1) + .2 * diff(yBounds);
                    else
                        whereY = yBounds(2) - .2 * diff(yBounds);
                    end                    
                    set(kidLines(length(kidLines) - 2), 'position', [xData(whereX) whereY], 'string', [' \bf' num2str(xData(whereX)) ', ' num2str(yData(whereX))]);
                end
            else
                for index = 1:length(kids) - 2
                    kidLines = get(kids(index), 'children');
                    set(kidLines(length(kidLines) - 1), 'xData', [xData(whereX) xData(whereX)], 'yData', get(kids(index), 'ylim'));
                    yData = get(kidLines(length(kidLines) - 4), 'yData');
                    firstPos = get(kidLines(length(kidLines) - 2), 'position');
                    yBounds = get(kids(index), 'ylim');
                    if abs(yBounds(1) - yData(whereX)) > abs(yBounds(2) - yData(whereX))
                        whereY = yBounds(1) + .1 * diff(yBounds);
                    else
                        whereY = yBounds(2) - .1 * diff(yBounds);
                    end
                    set(kidLines(length(kidLines) - 3), 'position', [xData(whereX) whereY], 'string', [' \bf \Delta' num2str(xData(whereX)- firstPos(1)) ', \Delta' num2str(yData(whereX) - yData(1 + round(firstPos(1) / diff(xData(1:2)))))], 'visible', 'on');                    
                end
            end            
        elseif pointerLoc(1) <= imageLoc(1) && pointerLoc(2) > 0.05
            set(gcf, 'pointer', 'top');
        elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) <= .05
            set(gcf, 'pointer', 'right');
        end
        if pointerLoc(2) < .015
            set(gcf, 'pointer', 'arrow');
        end
    end
end
            
    
function mouseUpScope(varargin)
    pointerLoc = get(gcf, 'CurrentPoint');
    imageLoc = get(gca, 'Position');
    if length(get(gcf, 'userData')) > 1    
        info = get(gcf, 'userData');
        set(gcf, 'userData', 0);
        blueRect = get(gcf, 'children');
        delete(blueRect(1));
        kids = get(gcf, 'children');
        firstKid = get(kids(1), 'children');
        kidData = get(firstKid(end - 4), 'xData');
        if info(1) == -1
            if info(2) ~= pointerLoc(1)
                myPoint = (pointerLoc(1) - .05) / .95  * diff(get(kids(1),'Xlim')) + min(get(kids(1),'Xlim'));
                if myPoint < min(kidData)
                    myPoint = min(kidData);
                end
                if myPoint > max(kidData)
                    myPoint = max(kidData);
                end                
                set(kids(1:length(kids) - 2), 'xlim', sort([(info(2) - .05) / .95  * diff(get(kids(1),'Xlim')) + min(get(kids(1),'Xlim')) myPoint]));
                set(kids(length(kids) - 1), 'string', num2str((max(kidData) - min(kidData)) / diff(get(kids(1), 'Xlim')),'%-1.1f'));
                xBounds = get(kids(1), 'xlim');
                
                zoomFactor = str2double(get(kids(length(kids) - 1), 'string'));
                newStep = 1 / zoomFactor / (1- 1 / zoomFactor);
                if newStep > 10
                    set(kids(length(kids)), 'sliderStep', [1 newStep]);
                else
                    set(kids(length(kids)), 'sliderStep', [newStep / 10 newStep]);
                end
                set(kids(length(kids)), 'value', xBounds(1) / max(kidData) / (1- 1 / zoomFactor));
            end
        else
            if info(2) ~= pointerLoc(2)
                set(kids(length(kids) - 1 - info(1)), 'ylim', sort((([info(2) pointerLoc(2)] - .05 - (info(1) - 1) * (.95 / (length(kids) - 2))) / (.95 / (length(kids) - 2)))  * diff(get(kids(length(kids) - 1 - info(1)),'ylim')) + min(get(kids(length(kids) - 1 - info(1)),'ylim'))));
            end
        end
    elseif strcmp(get(gcf, 'SelectionType'), 'alt')
        kids = get(gcf, 'children');
        if pointerLoc(1) <= imageLoc(1) && pointerLoc(2) > 0.05
            whichAxis = fix((pointerLoc(2) - .05) / .95 * (length(kids) - 2)) + 1;
            lineKids = get(kids(length(kids) - 1 - whichAxis), 'children');
            tempXData = get(lineKids(1), 'xData');
            tempYData = get(lineKids(1), 'yData');
            tempBounds = get(kids(length(kids) - 1 - whichAxis), 'xlim');
            minBound = find(tempXData == min(tempBounds));
            maxBound = find(tempXData == max(tempBounds));
            set(lineKids(length(lineKids)), 'yData', [min(tempYData(minBound:maxBound)) max(tempYData(minBound:maxBound))]);
            set(lineKids(length(lineKids) - 1), 'yData', [min(tempYData(minBound:maxBound)) max(tempYData(minBound:maxBound))]);
            set(kids(length(kids) - 1 - whichAxis), 'YLimMode', 'auto');
        elseif pointerLoc(1) > imageLoc(1) && pointerLoc(1) < imageLoc(1) + imageLoc(3) && pointerLoc(2) <= .05
            set(kids(1:length(kids) - 2), 'XLimMode', 'auto');
            set(kids(length(kids) - 1), 'string', '1.0');
            set(kids(length(kids)), 'sliderStep', [1 1/0]);
        end      
    end
end
    
    
function traceScroll(varargin)

    kids = get(gcf, 'children');
    firstKid = get(kids(1), 'children');
    xData = get(firstKid(end - 4), 'xdata');
    zoomFactor = str2double(get(kids(length(kids) - 1), 'string'));
    scrollValue = get(kids(length(kids)), 'value') * (1- 1 / zoomFactor);
    
    windowSize = (xData(end) - xData(1)) / zoomFactor;
    newStep = 1 / zoomFactor / (1- 1 / zoomFactor);
    if newStep > 10
        set(kids(length(kids)), 'sliderStep', [1 newStep]);
    else
        set(kids(length(kids)), 'sliderStep', [newStep / 10 newStep]);
    end
    set(kids(1:length(kids) - 2), 'Xlim', [xData(end) * scrollValue xData(end) * scrollValue + windowSize], 'dataaspectratiomode', 'auto', 'plotboxaspectratiomode', 'auto');
end
