%{
->simf.TraceParams      
->simf.RFResponses
---
trace                             : mediumblob                            #
%}


classdef Traces < dj.Computed
    properties
         keySource  = simf.TraceParams * simf.RFResponses
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            trace = fetch1(simf.RFResponses & key,'response');
            activ_func = getActivation(simf.TraceParams & key);
            key.trace = activ_func(trace);
            self.insert(key)
        end
    end
end

