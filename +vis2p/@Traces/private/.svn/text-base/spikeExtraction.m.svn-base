function [traces,out] = spikeExtraction(traces,fps,traceOpt)

traces = traces(:,:,1);
traces = getCaEvents( traces, fps, traceOpt.tau, traceOpt.highPass );


if nargout>1
    out = [];
end