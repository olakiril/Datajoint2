function [traces, qual, traceOpt] = oopsiBest(traces,fps,traceOpt)

traces = traces(:,:,1);
hp = traceOpt.highPass;
traces = traces + abs(min(traces(:)))+eps;
traces = traces./convmirr(traces,hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
traces = bsxfun(@plus,traces,abs(min(traces)))+eps;

if fps>20
    otraces = trresize(traces,fps,50,'linear');
    fps = 20;
    reduce = true;
else
    reduce = false;
    otraces = traces;
end

qual = nan(size(otraces,2),1);

ops = struct('lam',traceOpt.lamda);
if  ~isnan(traceOpt.tau)
    ops.gam = (1-1/fps/traceOpt.tau);
end
    
for iTrace = 1:size(traces,2)
    if  ~isnan(traceOpt.sigma)
        ops.sig = mean(mad( traces(:,iTrace)',1)*1.4826*traceOpt.sigma);
    end
    t = otraces(:,iTrace);
    [otraces(:,iTrace),~,~,d]= fast_oopsi(otraces(:,iTrace)', struct('dt',1/fps),ops);
    qual(iTrace) = corr(t,d);
end

if reduce
    xi = 0.5:size(otraces,1)/size(traces,1):size(otraces,1)+0.4;
    traces = interp1(1:size(otraces,1),otraces,xi);
else
    traces = otraces;
end

%% pre-2013-03-14
% function traces = oopsiBest(traces,fps,traceOpt)
% 
% traces = traces(:,:,1);
% hp = traceOpt.highPass;
% traces = traces + abs(min(traces(:)))+eps;
% traces = traces./convmirr(traces,hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
% traces = bsxfun(@plus,traces,abs(min(traces)))+eps;
% 
% otraces = trresize(traces,fps,100);
% fps = 10;
% for iTrace = 1:size(traces,2)
%     otraces(:,iTrace)= fast_oopsi(otraces(:,iTrace)', struct('dt',1/fps),...
%         struct('lam',traceOpt.lamda,'gam',(1-1/fps/traceOpt.tau)));
% end
% 
% xi = 0.5:size(otraces,1)/size(traces,1):size(otraces,1)+0.4;
% traces = interp1(1:size(otraces,1),otraces,xi);

