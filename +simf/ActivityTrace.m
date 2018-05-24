%{
->simf.Activity
filter_id                             : smallint      # 
---
trace                             : mediumblob                            #
%}


classdef ActivityTrace < dj.Part
    properties(SetAccess=protected)
        master = simf.Activity
    end

    methods
        function makeTuples(self, key)
            self.insert(key)
        end
    end
end

