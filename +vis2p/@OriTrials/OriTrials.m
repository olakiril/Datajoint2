%{
vis2p.OriTrials (computed) # 
-> vis2p.VisStims
repeat_num      : mediumint unsigned     # 
ori_num         : mediumint unsigned     # i) the orientation shown in this trial
---
ori_times                   : mediumblob                    # c) frame timestamps in win time 
%}


classdef OriTrials < dj.Relvar & dj.AutoPopulate

% 	properties
% 		popRel = vis2p.VisStims('exp_type = "GratingExperiment"  or exp_type = "MultDimExperiment"')
% 	end

	methods(Access=protected)

		makeTuples( obj, key)
    end
    
    methods
		function self = OriTrials(varargin)
			self.restrict(varargin{:})
		end


	end

end