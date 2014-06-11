%{
vis2p.StatsModel (computed) # 
-> vis2p.StatsSites
class_opt       : tinyint unsigned       # classification option
---
performance                 : mediumblob                    # c) distortion of the fit of the true data to the k-clusters
class_mean                  : mediumblob                    # c) mean of the responses
class_var                   : mediumblob                    # c) variance of the responses
class_per                   : mediumblob                    # c) distortion of the fit of the true data to the k-clusters
INDEX(mouse_id,exp_date,scan_idx,movie_type)
INDEX(mouse_id,exp_date,scan_idx,trace_opt,movie_type,stats_opt,stim_idx)
%}


classdef StatsModel < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.StatsSites('trace_opt = 17')*vis2p.StatsModelParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
    end
    
    methods
		function self = StatsModel(varargin)
			self.restrict(varargin{:})
		end

		Trace = getTraces( obj )


		plot(obj,key,varargin)

		[J, rJ] = run( obj, classes,classreps )

	end

end