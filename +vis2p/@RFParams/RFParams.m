%{
vis2p.RFParams (computed) #
-> vis2p.VisStims
dot_size        : smallint unsigned      # i) size of stimulus in pixels
stim_frames     : smallint unsigned      # i) length of stimulus in frames
---
loc_total_num               : smallint unsigned             # i) number of different locations
total_repetitions           : smallint unsigned             # i) total number of repetitions
%}


classdef RFParams < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        popRel = vis2p.VisStims('exp_type = "DotMappingExperiment" and process = "yes" or exp_type = "MouseDotMapping" and process = "yes"')
    end
    
    methods(Access=protected)
        
        makeTuples( obj, key )
            
    end
    
    methods
        function self = RFParams(varargin)
            self.restrict(varargin{:})
        end
        
        
    end
    
end