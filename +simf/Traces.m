%{
->simf.TraceGroup     
---
->simf.RFResponses
trace                             : mediumblob                            #
%}


classdef Traces < dj.Computed
   
    methods(Access=protected)
        function makeTuples(self, key)
            self.insert(key)
        end
    end
end

