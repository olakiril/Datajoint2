%{
vis2p.StatsAwake (computed) # 
-> vis2p.TracesGroup
-> vis2p.StatsAwakeParams
movie_type      : varchar(10)            # the type of movie shown
stats_opt       : smallint               # 
state           : varchar(10)            # 
---
mn                          : float(10,5)                   # c)average eucleadian distance
vr                          : float(10,5)                   # 
ps                          : float(10,5)                   # 
ls                          : float(10,5)                   # c) variance
cr                          : float(10,5)                   # c) population sparseness
rl                          : float(10,5)                   # c) probability of zero neurons responding
av_trials                   : float(10,5)                   # c) kurtosis
mi                          : float(10,5)                   # c) lifetime sparseness
INDEX(mouse_id,exp_date,scan_idx,movie_type)
INDEX(mouse_id,exp_date,scan_idx,stats_opt,trace_opt,movie_type)
INDEX(movie_type)
INDEX(stats_opt)
%}


classdef StatsAwake < dj.Relvar & dj.AutoPopulate

	properties
		popRel = (vis2p.Scans(' exp_date>"2013-09-17" and state= "awake" and problem_type = "none!"').*vis2p.StatsSites)*vis2p.StatsAwakeParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsAwake(varargin)
			self.restrict(varargin{:})
		end


	end

end