%{
vis2p.StatsSiteTracesParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
rf_snr_thr      : double(5,0) unsigned   # m) rf snr threshold
trace_qual      : double(5,0) unsigned   # m) trace quality threshold from skewness of the trace
rf_p_thr        : double(5,0) unsigned   # m) rf p threshold
---
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsSiteTracesParams < dj.Relvar
	methods

		function self = StatsSiteTracesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end