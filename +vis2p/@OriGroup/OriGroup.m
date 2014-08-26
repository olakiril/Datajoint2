%{
vis2p.OriGroup (computed) # 
-> vis2p.VisStims
---
uni_ori                     : mediumint unsigned            # 
total_trials                : mediumint unsigned            # i) the orientation shown in this trial
%}


classdef OriGroup < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.VisStims('exp_type = "GratingExperiment"  or exp_type = "MultDimExperiment"')
	end

	methods(Access=protected)

		makeTuples( obj, key)
    end
    
    methods
		function self = OriGroup(varargin)
			self.restrict(varargin{:})
		end


	end

end
