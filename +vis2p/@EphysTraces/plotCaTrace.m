function out = plotCaTrace(obj,varargin)

params.trace_opt = [3 20];
params.colors = [0.2 0.2 1;0.2 1 0.2;1 0.2 0.2];
params.names = [];

params = getParams(params,varargin);

import vis2p.*
if isempty(params.names)
    for i = 1:length(params.trace_opt)
        params.names{i} = num2str(params.trace_opt(i));
    end
end
params.names{end+1} = 'Spikes';

keys = fetch(obj);

if length(params.colors)<length(keys)
    params.colors = hsv(length(params.trace_opt) +1);
end

fh = figure;
set(fh,'position',[626 319 1042 592])
for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
%     if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        spikeTimes = fetch1(EphysTraces(key),'spike_times');
        Fs = fetch1(Movies(key),'fps');
      
        key.masknum = fetchn(EphysTraces(key),'masknum');
       
        hold on
        
        for itraceopt = 1:length(params.trace_opt);
            key.trace_opt = params.trace_opt(itraceopt);
            trace = fetch1(Traces(key),'trace');
            if params.trace_opt(itraceopt) == 3
                trace = trace;
            end
            plot(0:1/Fs:(length(trace)-1)/Fs,trace,'color',params.colors(itraceopt,:));
            
        end
        
        plot(spikeTimes/1000,zeros(length(spikeTimes),1),'.','Color',params.colors(end,:))

        title([key.exp_date ' scan:' num2str(key.scan_idx)])
        xlabel('Time (sec)')
        ylabel('DF/F')
        l = legend(params.names);
        legend(l,'BoxOff');
        set(fh,'name',[key.exp_date ' ' num2str(key.scan_idx)])
        
        if ikey ~= 1
            pause
        end
%     end
    if nargout>0
        out = fh;
    end
end