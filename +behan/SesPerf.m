%{
# Performance with 2AFC tasks
->beh.Session
---
%}

classdef SesPerf < dj.Computed
    properties
        keySource = beh.Session & (beh.RewardCond & 'probe=1') & (beh.RewardCond & 'probe=2') & 'animal_id>0'
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            insert(self,key)
            
        end
    end
    
    methods
        function updateTrials(self,restrict)
            
            % update session key
            populate(self,restrict)
            
            % update trials
            populate(behan.Perf,restrict)
            
        end
    end
end