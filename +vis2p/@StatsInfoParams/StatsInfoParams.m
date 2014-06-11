%{
vis2p.StatsInfoParams (lookup) # 
binsize         : smallint unsigned      # m) ms binsize for correlation computations
neurons         : smallint unsigned      # m) number of neurons in each cluster
thr_factor      : smallint unsigned      # m) factor to multiply the std to find  threshold 
edge_crop       : smallint unsigned      # m) size of trace chopped from each trial (ms)
dist_fact       : enum('3','0.3','0','-0.3','-3','random')  # m) cluster selection: factorization of distances
trace_opt       : smallint               # 
---
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsInfoParams < dj.Relvar
	methods

		function self = StatsInfoParams(varargin)
			self.restrict(varargin{:})
		end
	end

end