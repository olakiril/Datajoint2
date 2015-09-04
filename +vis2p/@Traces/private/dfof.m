function [traces,out] = dfof(traces,fps,traceOpt)

traces = traces(:,:,1);
traces = traces+495424;
if isfield(traceOpt,'highPass') && ~isnan(traceOpt.highPass)
    if traceOpt.highPass>0
        k = hamming(round(fps/traceOpt.highPass)*2+1);
        k = k/sum(k);
        traces = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass
    else
        traces = bsxfun(@rdivide,traces,mean(traces))-1;  % dF/F where F is mean
    end
end
% low-pass filtration
if isfield(traceOpt,'lowPass') && ~isnan(traceOpt.lowPass)
    k = hamming(round(fps/traceOpt.lowPass)*2+1);
    k = k/sum(k);
    traces = convmirr(traces,k);
end


if nargout>1
    out = zeros(size(traces,1));
end