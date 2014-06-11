%{
vis2p.CenterSurTrials (computed) #
mouse_id        : smallint unsigned      #
exp_date        : date                   #
scan_idx        : smallint unsigned      #
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
repeat_num      : mediumint unsigned     #
---
ori_in                      : mediumint                     # i) the number of the movie shown in this trial
ori_out                     : mediumint                     # i) the number of the movie shown in this trial
movie_times                 : mediumblob                    # c) movieframe timestamps in win time
%}


classdef CenterSurTrials < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        popRel = vis2p.VisStims('exp_type = "CenterSurround"')
    end
    
    methods
        function self = CenterSurTrials(varargin)
            self.restrict(varargin{:})
        end
        
    end
    methods(Access=protected)
        
        
        makeTuples( obj, key)
        
    end
    
end