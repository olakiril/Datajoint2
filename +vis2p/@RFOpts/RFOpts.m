%{
vis2p.RFOpts (lookup) # 
rf_opt_num      : smallint unsigned      # option index
---
time_on                     : mediumint                     # m) time that the response window starts relative to the stimulus (ms)
time_off                    : mediumint                     # m) time that the response window ends relative to the stimulus (ms)
fraction_of_data            : float(6,2) unsigned           # m) the fraction of data that should be used for the analysis
%}


classdef RFOpts < dj.Relvar
	methods

		function self = RFOpts(varargin)
			self.restrict(varargin{:})
		end
	end

end