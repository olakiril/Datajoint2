%{
vis2p.StatsRf (computed) # 
-> vis2p.TracesGroup
-> vis2p.StatsRfParams
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
movie_type      : varchar(10)            # the type of movie shown
---
actual_binsize              : float                         # c) (ms) actual binsize used
pzero_site                  : mediumblob                    # c) probability of zero neurons responding
%}


classdef StatsRf < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.TracesGroup('trace_opt = 6')*vis2p.StatsParams('movie_type = "phase"')*vis2p.StatsRfParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsRf(varargin)
			self.restrict(varargin{:})
		end

		data = getTraceStruct(obj, key )

		traces = getTraces(obj, key )


	end

end