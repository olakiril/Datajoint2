%{
vis2p.StatsSitesParams (lookup) # 
stats_opt       : smallint               # 
---
binsize=100                 : smallint unsigned             # m) ms binsize for correlation computations
rf_snr_thr=null             : double(5,0) unsigned          # m) rf snr threshold
trace_qual=null             : double(5,3) unsigned          # m) trace quality threshold from skewness of the trace
rf_p_thr=null               : double(5,0) unsigned          # m) rf p threshold
event_thr=null              : double(3,0)                   # m) % of std of noise for event characterization
samp_method=null            : varchar(10)                   # 
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsSitesParams < dj.Relvar
	methods

		function self = StatsSitesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end