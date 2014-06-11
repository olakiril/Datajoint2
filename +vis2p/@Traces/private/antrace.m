function [traces,out] = antrace(traces,fps,traceOpts) %#ok<*INUSD>

traces = traces(:,:,3);


if nargout>1
    out = [];
end