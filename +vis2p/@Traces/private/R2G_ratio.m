function [traces, out, traceOpt] = R2G_ratio(traces,fps,traceOpt)

traces = traces(:,:,2)./traces(:,:,1);

if isfield(traceOpt,'highPass')
    if traceOpt.highPass>0
        k = hamming(round(fps/traceOpt.highPass)*2+1);
        k = k/sum(k);
        % make sure everything is positive
        traces = traces + abs(min(traces(:)))+eps;
        traces = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass
        
        % whatever is nan get rid of it by replacing it with 0
        traces(isnan(traces)) = 0;
    else
        traces = bsxfun(@rdivide,traces,mean(traces))-1;  % dF/F where F is mean
    end
else
    % make sure everything is positive
    traces = traces + abs(min(traces(:)))+eps;
end

if nargout>1
    out = zeros(size(traces,1),1);
end