%{
vis2p.OriTraces (computed) # 
-> vis2p.Traces
-> vis2p.OriTracesParams
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
---
ori                         : mediumblob                    # c) the raw data [trials,oris]
Poti                        : double(4,3)                   # c) the significance of tuning from oti
Pdoti                       : double(4,3)                   # c) the significance of tuning from dprime
Pdm                         : mediumblob                    # c) oreffered direction of motion
INDEX(mouse_id,exp_date,scan_idx,stim_idx)
%}


classdef OriTraces < dj.Relvar & dj.AutoPopulate

	properties(Constant)
		popRel = Traces*OriTracesParams('process = "yes"')*VisStims('exp_type = "GratingExperiment" or exp_type = "MultDimExperiment"')
	end

	methods(Access=protected)
		makeTuples(self, key)
    end
    
    methods
		function self = OriTraces(varargin)
			self.restrict(varargin{:})
		end
	end

end