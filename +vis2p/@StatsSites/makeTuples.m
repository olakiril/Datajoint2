function makeTuples( obj, key )

% get traces
traces = getTraces(StatsSites,'key',key,'compute',1);

% control error
if isempty(traces)
    display('Nothing to compute!')
    return
end

mtraces = mean(traces,3);
ttraces = traces(:,:);

% trace info     
bin = fetch1(StatsSitesParams(key),'binsize');
key.time = size(traces,2)*bin/1000; % trace size in seconds;
key.neurons = size(traces,1);
key.trials = size(traces,3);

% mean
key.mean = nanmean(traces(:));

% Variance explained
key.varexp = reliability(traces);

% variance
key.variance = nanmean(reshape(var(traces,[],3),[],1));

% signal correlations
c = double(corr(mtraces'));
key.sigcorr = nanmean(c(logical(tril(ones(size(c)),-1))));
        
% pdist
key.eudist = nanmean(pdist(mtraces'));

% compute synchrony
v = (mtraces - mtraces.^2).^2;
v = mean(v,2);
m = mean(mtraces,2);
mv = (m - m.^2).^2;
key.synchrony = nanmean(mv./v);

% Population sparseness, adapted from Vinje & Galland 2000
key.pspars = nanmean(sparseness(mtraces));
key.pspars_tr =  nanmean(sparseness(ttraces));
key.lspars = nanmean(sparseness(mtraces'));
key.lspars_tr =  nanmean(sparseness(ttraces'));

% probability of zero response 
key.pzero = nanmean(sparseness(mtraces,'type','pzero'));
key.pzero_tr = nanmean(sparseness(ttraces,'type','pzero'));
key.lzero = nanmean(sparseness(mtraces','type','pzero'));
key.lzero_tr = nanmean(sparseness(ttraces','type','pzero'));

% kurtosis
key.pkurt = nanmean(fixinf(sparseness(mtraces,'type','kurtosis')));
key.pkurt_tr = nanmean(fixinf(sparseness(ttraces,'type','kurtosis')));
key.lkurt = nanmean(fixinf(sparseness(mtraces','type','kurtosis')));
key.lkurt_tr = nanmean(fixinf(sparseness(ttraces','type','kurtosis')));

% insert data
insert( obj, key );

function x=fixinf(x)
x(isinf(x)) = 0;

