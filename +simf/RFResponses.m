%{
->simf.RFRespGroup
->simf.RFFilters
---
response                             : mediumblob                            #
%}


classdef RFResponses < dj.Computed
    methods(Access=protected)
        function makeTuples(self, key)
                insert(self,key)
        end
    end
end

