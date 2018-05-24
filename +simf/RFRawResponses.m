%{
->simf.RFRaw
filter_id                             : smallint      # 
---
response                             : mediumblob                            #
%}


classdef RFRawResponses < dj.Part
    properties(SetAccess=protected)
        master = simf.RFRaw
    end
    methods
        function makeTuples(self, key)
                insert(self,key)
        end
    end
end

