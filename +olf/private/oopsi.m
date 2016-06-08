function [traces, qual, traceOpt] = oopsi(traces,fps,traceOpt)

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
 fitf = @(x) 0.4723*exp(-((x-16.46)/5.699)^2) + 0.3295*exp(-((x-9.802)/3.988)^2);
qual = nan(size(traces,2),1);
parfor iTrace = 1:size(traces,2)
    t = traces(:,iTrace);
    if  traceOpt.lamda<0
        k = kurtosis(traces(:,iTrace));
         ops = struct('lam',abs(traceOpt.lamda)+fitf(k)*1.75);
    else
        ops = struct('lam',traceOpt.lamda);
    end
    if  ~isnan(traceOpt.sigma)
        ops.sig = mean(mad( traces(:,iTrace)',1)*1.4826*traceOpt.sigma);
    end

    [traces(:,iTrace),~,~,d] = fast_oopsi( traces(:,iTrace)', struct('dt',1/fps),ops);
    qual(iTrace) = corr(t,d);
end

