function makeTuples( this, key )

% get traces
[cellnums,traces] = fetchn( vis2p.MaskTraces(key,'masknum>0'), 'masknum', 'calcium_trace' );
traces = [traces{:}];
fps    = fetch1( vis2p.Movies(key), 'fps' );
% filter traces
[etraces,cleanTraces,baseline] = getCaEvents(double(traces),fps); %#ok<ASGLU>

% compute snr
snr = std( cleanTraces ) ./ std( traces - cleanTraces - baseline );
for iCell=1:length(cellnums)
    tuple = key;
    tuple.masknum = cellnums(iCell);
    tuple.ca_event_snr = snr(iCell);
    insert( this, tuple );
end
