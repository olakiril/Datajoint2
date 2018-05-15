%{
->simf.TraceParams      
->simf.RFRespGroup
---
%}


classdef TraceGroup < dj.Computed
    properties
         keySource  = simf.TraceParams * simf.RFRespGroup
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            [traces,keys] = fetchn(simf.RFResponses & key,'response');
            traces = cell2mat(traces);
            activ_func = getActivation(simf.TraceParams & key);
            act_traces = activ_func(traces);
            [keys.trace_opt] = deal(key.trace_opt);
            
            insert(self,key)
            for ikey = 1:length(keys)
                trace_key = keys(ikey);
                trace_key.trace = act_traces(ikey,:);
                makeTuples(simf.Traces,trace_key)
            end
        end
    end
end

