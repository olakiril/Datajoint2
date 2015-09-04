%{
vis2p.OriResponse (computed) # 
-> vis2p.Traces
-> vis2p.OriTracesParams
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
---
Poti                        : double(4,3)                   # c) the significance of tuning from oti
Pdm                         : double                        # c) oreffered direction of motion
fitVM                       : mediumblob                    # Von misses fit
INDEX(mouse_id,exp_date,scan_idx,stim_idx)
%}


classdef OriResponse < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Traces('masknum>0') ... 
            *(vis2p.OriGroup('uni_ori<7').*vis2p.Scans('problem_type = "none!"'))
	end

	methods(Access=protected)
		makeTuples(self, key)
    end
    
    methods
		function self = OriResponse(varargin)
			self.restrict(varargin{:})
        end
        
        [T,oris] = getTraces(obj,varargin)
        
        plot(obj,type)
	end

end