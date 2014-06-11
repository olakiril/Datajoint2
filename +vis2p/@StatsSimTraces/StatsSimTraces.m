%{
vis2p.StatsSimTraces (computed) # 
movie_num       : varchar(20)            # the number of the movie shown
movie_type      : varchar(10)            # the type of movie shown
sim_traces_opt  : smallint unsigned      # 
---
sim_traces                  : mediumblob                    # c) traces from simulated RFs
%}


classdef StatsSimTraces < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.StatsSimTracesParams('process = "yes"')
	end

	methods(Access=protected)
        makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsSimTraces(varargin)
			self.restrict(varargin{:})
		end

		filters = getFilters(obj,key,varargin) 	


	end

end