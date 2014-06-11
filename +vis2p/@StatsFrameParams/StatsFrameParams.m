%{
vis2p.StatsFrameParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
delay           : smallint               # m) calcium response delay
---
process="yes"               : enum('no','yes')              # m) compute or not compute
INDEX(binsize)
%}


classdef StatsFrameParams < dj.Relvar
	methods

		function self = StatsFrameParams(varargin)
			self.restrict(varargin{:})
		end
	end

end