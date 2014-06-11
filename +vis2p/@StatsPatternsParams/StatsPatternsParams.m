%{
vis2p.StatsPatternsParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
thr_factor      : smallint unsigned      # m) factor to multiply the std to find  threshold 
---
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsPatternsParams < dj.Relvar
	methods

		function self = StatsPatternsParams(varargin)
			self.restrict(varargin{:})
		end
	end

end