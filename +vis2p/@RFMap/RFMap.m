%{
vis2p.RFMap (computed) # 
-> vis2p.Traces
-> vis2p.RFParams
-> vis2p.RFOpts
---
repeat_avg                  : smallint unsigned             # c) number of repeats that are averaged
on_rf=null                  : mediumblob                    # c) mean activity for each location
off_rf=null                 : mediumblob                    # c) normalized mean activity across all locations
onmoff_rf=null              : mediumblob                    # c) p value for each loaction
onpoff_rf=null              : mediumblob                    # c) p value for each loaction
INDEX(mouse_id,exp_date,scan_idx,dot_size,stim_frames)
INDEX(mouse_id,exp_date,scan_idx,masknum,dot_size,stim_frames,repeat_avg)
%}


classdef RFMap < dj.Relvar & dj.AutoPopulate

	properties(Constant)
		popRel = (vis2p.RFParams*vis2p.Traces)*vis2p.RFOpts
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = RFMap(varargin)
			self.restrict(varargin{:})
		end


	end

end