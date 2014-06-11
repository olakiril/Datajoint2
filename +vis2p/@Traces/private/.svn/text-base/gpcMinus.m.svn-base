function [traces,out,traceOpt] = gpcMinus(traces,fps,traceOpt)


[u,~] = eigs( cov(double(traces(:,:,3))),1 );
traces = double(traces(:,:,1));
traces = (traces - traces*u*u');

if nargout>1
    out = [];
end