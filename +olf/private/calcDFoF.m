function [traceDFoF,out,traceOpts] = calcDFoF(traces,Fps,traceOpts) %#ok<INUSD>

% function traceDFoF = calcDFoF(trace)
%
% calculates DFoF and does some minor baseline shift correction
%
% MF 2010-11-06

% get rid over any other traces
traces = traces(:,:,1);
traces = traces - min(traces(:));
params.WindowBaseline = 15;
params.Percentile = 10;
traceDFoF = nan(size(traces));
for itrace = 1:size(traces,2)
    trace = reshape(traces(:,itrace), 1, []);
    winSize = ceil(params.WindowBaseline*Fps);
    if length(trace)>winSize
        
        % Cut the trace into windows
        if mod(length(trace), winSize) ~= 0
            nbRepElements = winSize - mod(length(trace), winSize);
            traceWin = [trace repmat(trace(end), 1, nbRepElements)];
        else
            traceWin = trace;
        end
        
        traceWin = reshape(traceWin, winSize, []);
        segMin = prctile((traceWin)',10,2);

        % Prepare interpolation
        segPos = (winSize/2) + winSize * (0:length(segMin)-1);
        segPos(1) = 1;
        segPos(end) = length(trace);
        
        baseline = interp1(segPos, segMin, 1:length(trace), 'spline');
        
        Trace = (trace - (baseline-segMin(1)))';
    else
        Trace = trace';
        
    end
    
    % Calculate DF/F
    traceDFoF(:,itrace) = ((Trace - (prctile(Trace,params.Percentile)))/prctile(Trace,params.Percentile));
end


if nargout>1
    out = zeros(size(traces,1),1);
end