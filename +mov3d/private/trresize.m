function [trace, d] = trresize(trace,fps,bin,method)



if nargin<4 || strcmp(method,'conv')
    d = max(1,round(bin/1000*fps));
    k = ones(d,1)/d;
    trace = convn(trace,k,'same');
    trace = trace(1:d:end,:,:);
elseif strcmp(method,'bin')
    d = max(1,round(bin/1000*fps));
    trace(isnan(trace)) = 0;
    trace = cumsum(trace);
    trace = trace(1:d:end,:,:);
    trace = (trace(2:end,:,:) - trace(1:end-1,:,:))/d;
elseif strcmp(method,'binsum')
    d = max(1,round(bin/1000*fps));
    trace(isnan(trace)) = 0;
    trace = cumsum(trace);
    trace = trace(1:d:end,:,:);
    trace = (trace(2:end,:,:) - trace(1:end-1,:,:));
elseif strcmp(method,'hamming')
    d = bin/1000*fps;
    trace = interp1(trace,d/4:d/4:size(trace,1));
    if size(trace,1)==1
        trace = trace';
    end
    d = 4;
    k = hamming(d*2+1);
    k = k/sum(k);
    trace = single(convmirr(double(trace),k));
    trace = trace(ceil(d/2):d:end,:,:);
else
    d = bin/1000*fps;
    trace = interp1(1:size(trace,1),trace,1:d:size(trace,1),method);
    if size(trace,1)==1;trace = trace';end
end
