function makeTuples( obj, key )

areas = {'v1i' 'v1o' 'v2i' 'v2o'};

for iarea = 1:length(areas)
    display(areas{iarea})
    tuple = key;
    tuple.area = areas{iarea};
    
    % get traces
    traces = getTraces(StatArea(key),'area',areas{iarea},'key',key);
    
    % get masknums
    maskIds = fetch1(StatArea(key),areas{iarea});
    maskKey = sprintf('masknum=%d or ',maskIds);
    maskKey = maskKey(1:end-3);
    
    % control error
    if isempty(traces)
        display('Nothing to compute!')
        return
    end
    
    mtraces = mean(traces,3);
    
    % trace info
    bin = fetch1(StatsSitesParams(tuple),'binsize');
    tuple.time = size(traces,2)*bin/1000; % trace size in seconds;
    tuple.neurons = size(traces,1);
    tuple.trials = size(traces,3);
    
    % mean
    tuple.mean = nanmean(traces(:));
    
    % variance
    tuple.variance = nanmean(reshape(var(traces,[],3),[],1));
    
    % signal correlations
    c = double(corr(mtraces'));
    tuple.sigcorr = nanmean(c(logical(tril(ones(size(c)),-1))));
    
    % Population sparseness, adapted from Vinje & Galland 2000
    tuple.pspars = nanmean(sparseness(mtraces));
    tuple.lspars = nanmean(sparseness(mtraces'));
    
    % classification performance
    tuple.pMulti= nnclassRaw(traces);
    tuple.pDist= nnclass(traces,'cells',size(traces,1));
    
    % reliability
    tuple.reliability = double(mean((fetchn(StatsTraces(key,maskKey),'inCorr'))));
    
    % insert data
    insert( obj, tuple );
    
end
