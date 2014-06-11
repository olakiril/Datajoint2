%{
vis2p.RFMapRaw (computed) # 
-> vis2p.RFParams
-> vis2p.RFOpts
---
onpoff_rf                   : longblob                      # c) p value for each loaction
INDEX(mouse_id,exp_date,scan_idx,dot_size,stim_frames)
INDEX(mouse_id,exp_date,scan_idx,dot_size,stim_frames)
INDEX(mouse_id,exp_date,scan_idx)
%}


classdef RFMapRaw < dj.Relvar & dj.AutoPopulate

	properties
		popRel = (vis2p.RFParams*vis2p.Scans)*vis2p.RFOpts('rf_opt_num = 3')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = RFMapRaw(varargin)
			self.restrict(varargin{:})
		end


		plot(obj,key)

	end

end