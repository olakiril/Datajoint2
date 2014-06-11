%{
vis2p.StatsTracesParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
undersample     : tinyint                # m) undersample after bining
event_thr       : tinyint                # m) the multiplication factor of std for the threshold calculation
---
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsTracesParams < dj.Relvar
	methods

		function self = StatsTracesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end