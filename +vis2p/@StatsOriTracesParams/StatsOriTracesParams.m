%{
vis2p.StatsOriTracesParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
undersample     : tinyint                # m) undersample after bining
shuffle         : mediumint              # m) number of shuffles
resp_delay      : smallint               # m) the dalay between stimulus and response (msec
---
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsOriTracesParams < dj.Relvar
	methods

		function self = StatsOriTracesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end