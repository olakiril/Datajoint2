%{
vis2p.StatsMulti (computed) # 
-> vis2p.StatsSites
multi_opt       : tinyint unsigned       # classification option
---
performance                 : double                        # c) classification performance
INDEX(mouse_id,exp_date,scan_idx,movie_type)
INDEX(mouse_id,exp_date,scan_idx,trace_opt,movie_type,stats_opt,stim_idx)
%}


classdef StatsMulti < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.TracesGroup*vis2p.StatsSites*vis2p.StatsMultiParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsMulti(varargin)
			self.restrict(varargin{:})
		end


	end

end