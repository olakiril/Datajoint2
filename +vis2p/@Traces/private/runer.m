keys = fetch(EphysTraces.*Scans('scan_prog = "MpScan"'));

parfor ikey = 1:length(keys)
    key = keys(ikey);
    key.trace_opt =17;
    populate(TracesGroup2,key);
end


%%
keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

parfor ikey = 1:length(keys)
    key = keys(ikey);
    key.trace_opt =27;
    populate(TracesGroup,key);
end

%%
keys = fetch(StatArea.*StatsSites('exp_date>"2013-08-05"'));

parfor ikey = 1:length(keys)
    key = keys(ikey);
    key.trace_opt =28;
%     populate(TracesGroup,key);
    populate(StatsSites,key);
    populate(StatArea,key);
    populate(StatAreaData,key);
end