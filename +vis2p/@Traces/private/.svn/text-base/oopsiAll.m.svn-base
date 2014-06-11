function [traces, qual,traceOpt] = oopsiAll(traces,fps,traceOpt)

traces = traces(:,:,1);

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

[traceAll,~,~,d] = fast_oopsi( traces(:)', struct('dt',1/fps),struct('lam',traceOpt.lamda));
tracesP = reshape(traceAll,size(traces));
tracesD = reshape(d,size(traces));

qual = nan(size(traces,2),1);
for iTrace = 1:size(traces,2)
    qual(iTrace) = corr(traces(:,iTrace),tracesD(:,iTrace));
end

traces = tracesP;

