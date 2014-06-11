%{
vis2p.StatsSim (computed) # 
movie_num       : varchar(20)            # the number of the movie shown
movie_type      : varchar(10)            # the type of movie shown
sim_traces_opt  : smallint unsigned      # 
sim_opt         : smallint unsigned      # 
---
p_spars                     : mediumblob                    # 
l_spars                     : mediumblob                    # 
mean                        : mediumblob                    # 
variance                    : mediumblob                    # 
e_dist                      : mediumblob                    # 
d_perf                      : mediumblob                    # 
s_corr                      : mediumblob                    # 
%}


classdef StatsSim < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.StatsSimTraces*vis2p.StatsSimParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsSim(varargin)
			self.restrict(varargin{:})
		end


	end

end