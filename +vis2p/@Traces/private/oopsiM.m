function [traces, qual] = oopsiM(traces,fps,traceOpt)

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

% Convolution
GaussWin = ceil(3*fps);
for iTrace = 1:size(traces,2)
    ctraces = conv(traces(:,iTrace) ,gausswin(GaussWin));   %convolution with the gaussian window
    ctraces = ctraces(floor(GaussWin/2):end-(ceil(GaussWin/2)));  %Correct for the shift
    traces(:,iTrace) = ctraces*prctile(traces(:,iTrace),80)/prctile(ctraces,80);
end

qual = nan(size(traces,2),1);
for iTrace = 1:size(traces,2)
    t = traces(:,iTrace);
    [traces(:,iTrace) ,~,~,d] = fast_oopsi( traces(:,iTrace)', struct('dt',1/fps));
    qual(iTrace) = corr(t,d);
end
for iTrace = 1:size(traces,2)
    traces(traces(:,iTrace)<3*std(traces(:,iTrace)),iTrace) = 0;
end
