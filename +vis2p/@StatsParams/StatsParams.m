%{
vis2p.StatsParams (computed) # 
-> vis2p.VisStims
movie_type      : varchar(10)            # the type of movie shown
trial_trigger   : tinyint                # 
---
INDEX(movie_type)
INDEX(movie_type,stim_idx)
%}


classdef StatsParams < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.VisStims('exp_type = "MouseStatExperiment" or exp_type = "MoviesExperiment"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsParams(varargin)
			self.restrict(varargin{:})
		end


	end

end