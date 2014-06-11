%{
vis2p.StatAreaDecodeRes (computed) # 
-> vis2p.StatAreaDecode
-> vis2p.StatAreaDecodeResParams
---
mutinfo                     : mediumblob                    # 
cperf                       : mediumblob                    # 
confmat                     : mediumblob                    # 
INDEX(dec_opt)
INDEX(stats_opt)
%}


classdef StatAreaDecodeRes < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.StatAreaDecode*vis2p.StatAreaDecodeResParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
    end
    
    methods
		function self = StatAreaDecodeRes(varargin)
			self.restrict(varargin{:})
		end


	end

end