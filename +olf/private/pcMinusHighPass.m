function [traces,out,traceOpt] = pcMinusHighPass(traces,fps,traceOpt) 

traces = double(traces(:,:,1));

if isfield(traceOpt,'highPass')
    if traceOpt.highPass>0
        
        % make sure everything is positive
        traces = traces + abs(min(traces(:)))+eps;
        
        % heavy filtering for pca removal
        k = hamming(round(fps/(traceOpt.highPass*20))*2+1);k = k/sum(k);
        tracesF = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass

        % normal filtering 
        k = hamming(round(fps/traceOpt.highPass)*2+1);k = k/sum(k);
        traces = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass
        
        % whatever is nan get rid of it by replacing it with 0
        traces(isnan(traces)) = 0; tracesF(isnan(tracesF)) = 0;

        % disable further filtering
        traceOpt = rmfield(traceOpt,'highPass');
    else
        traces = bsxfun(@rdivide,traces,mean(traces))-1;  % dF/F where F is mean
    end
end


[u,~] = eigs( cov(double(tracesF)),1 );
traces = (traces - traces*u*u');


if nargout>1
    out = [];
end


