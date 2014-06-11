%{
vis2p.MovieParams (computed) #
-> vis2p.VisStims
---
%}


classdef MovieParams < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = vis2p.VisStims('exp_type = "MouseMovieExperiment"')
    end
    
    methods(Access=protected)
        
        
        
        makeTuples( obj, key )
        
    end
    methods
        function self = MovieParams(varargin)
            self.restrict(varargin{:})
        end
    end
    
end