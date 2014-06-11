function [traces,out] = rtrace(traces,fps,traceOpts) %#ok<*INUSD>

traces = traces(:,:,2);

if nargout>1
    out = [];
end