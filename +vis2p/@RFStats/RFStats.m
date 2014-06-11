%{
vis2p.RFStats (computed) # 
-> vis2p.Traces
-> vis2p.RFFit
-> vis2p.RFParams
-> vis2p.RFOpts
---
onpoff_p                    : float                         # c) p value for each loaction
INDEX(mouse_id,exp_date,scan_idx,dot_size,stim_frames)
INDEX(mouse_id,exp_date,scan_idx,masknum,dot_size,stim_frames)
%}


classdef RFStats < dj.Relvar & dj.AutoPopulate

	properties(Constant)
		popRel = vis2p.RFFit
	end

	methods(Access=protected)

        
		makeTuples( obj, key )
        
    end
    
    methods
		function self = RFStats(varargin)
			self.restrict(varargin{:})
		end


	end

end