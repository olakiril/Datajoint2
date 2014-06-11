function makeTuples( obj, key )

classifier = fetch1(StatsMultiParams(key),'classifier');
traces = getTraces(StatsSites(key));

if isempty(traces)
    display('Nothing to compute!')
    return
end

eval(['P = ' classifier '(traces);']);

key.performance = P;

% insert data
insert( obj, key );
