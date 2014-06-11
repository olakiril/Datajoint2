function [traces,out] = traceFilter(traces,fps,traceOpt)

% high pass filtration
if isfield(traceOpt,'highPass') && ~isnan(traceOpt.highPass) && traceOpt.highPass>0
    k = hamming(round(fps/traceOpt.highPass)*2+1);
    k = k/sum(k);
    traces = convmirr(traces,k);  %  dF/F where F is low pass
end

% low-pass filtration
if isfield(traceOpt,'lowPass') && ~isnan(traceOpt.lowPass)
    k = hamming(round(fps/traceOpt.lowPass)*2+1);
    k = k/sum(k);
    traces = convmirr(traces,k);
end


if nargout>1
    out = [];
end