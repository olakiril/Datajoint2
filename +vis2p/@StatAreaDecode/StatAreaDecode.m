%{
vis2p.StatAreaDecode (computed) # 
-> vis2p.StatsSitesParams
-> vis2p.StatAreaDecodeParams
mouse_id        : smallint unsigned      # 
exp_date        : date                   # 
scan_idx        : smallint unsigned      # 
trace_opt       : smallint unsigned      # 
movie_type      : varchar(10)            # the type of movie shown
area_opt        : smallint               # 
---
cpV1                        : mediumblob                    # 
cpV2                        : mediumblob                    # 
%}


classdef StatAreaDecode < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.StatsSites('stats_opt = 7 and trace_opt = 17')*vis2p.StatArea*vis2p.StatAreaDecodeParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
    end
    
    methods
		function self = StatAreaDecode(varargin)
			self.restrict(varargin{:})
		end


	end

end