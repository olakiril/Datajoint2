function [regOUT pOUT rOUT] = regressCa(obj,varargin)

params.trace_opt = 17;
params.bin = 1000;
params.name = [];
params.plot = 1;

params = getParams(params,varargin);

if isempty(params.name)
    params.name = ['trace_opt: ' num2str(params.trace_opt)];
end

import vis2p.*
if nargout>0
    params.plot = 0;
end

keys = fetch(obj);

for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
    key.trace_opt = params.trace_opt;
    
    if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        
        % get traces
        fps = fetchn(Movies.*EphysTraces(key),'fps');
        
        traceP = fetchn(Traces(key,'masknum = 0'),'trace');
        key.masknum = fetchn(EphysTraces(key),'masknum');
        traceC  = getCaTraces(EphysTraces(key));
        
        traceC = trresize(traceC{1},fps,params.bin,'bin');
        traceP = trresize(traceP{1},fps,params.bin,'binsum')/params.bin*1000;
        
        [reg p r] = regressPlot(traceP,traceC,'thr',0.001,'plot',params.plot);
        
        if nargout>0
            regOUT = reg;
            pOUT = p;
            rOUT = r;
        else
            % set figure properties
            title([key.exp_date ' scan:' num2str(key.scan_idx)])
            xlabel('Firing Rate (sp/sec)')
            ylabel(params.name)
            set(gcf,'name',[key.exp_date ' ' num2str(key.scan_idx)])
        end
        
        if ikey ~= 1
            pause
        end
    end
end