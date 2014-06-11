
function out = plot(keys)


import vis2p.*
keys = fetch(keys);
fh = figure;
for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
    if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        eData = fetch1(MaskTraces(key,'masknum = 0'),'ephys_trace');
        spikeTimes = fetch1(EphysTraces(key),'spike_times');
        Fs = fetch1(Movies(key),'ephys_fs');
          
        spikeIdx = spikeTimes/1000*Fs + 1; % convert to index
        plot(eData);
        hold on;
        plot(spikeIdx,zeros(size(spikeIdx)),'.r');
        title([key.exp_date ' ' num2str(key.scan_idx)])
        
        if ikey ~= 1
            pause
        end
    end
    if nargout>0
        out = fh;
    end
end