function [traces, qual] = oopsiParam(traces,fps,traceOpt)

traces = traces(:,:,1);

if isfield(traceOpt,'highPass')
    if traceOpt.highPass>0
        k = hamming(round(fps/traceOpt.highPass)*2+1);
        k = k/sum(k);
        traces = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass
    else
        traces = bsxfun(@rdivide,traces,mean(traces))-1;  % dF/F where F is mean
    end
end
qual = nan(size(traces,2),1);
for iTrace = 1:size(traces,2)
    t = traces(:,iTrace);
    [traces(:,iTrace),~,~,d] = fast_oopsi( traces(:,iTrace)', struct('dt',1/fps),struct('lam',traceOpt.avgfr));
    qual(iTrace) = corr(t,d);
end

