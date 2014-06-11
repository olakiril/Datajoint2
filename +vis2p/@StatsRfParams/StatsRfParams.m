%{
vis2p.StatsRfParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
delay           : smallint               # m) calcium response delay
---
process="yes"               : enum('no','yes')              # m) compute or not compute
INDEX(binsize)
%}


classdef StatsRfParams < dj.Relvar
	methods

		function self = StatsRfParams(varargin)
			self.restrict(varargin{:})
		end
	end

end