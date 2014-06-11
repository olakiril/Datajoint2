%{
vis2p.StatsSimParams (lookup) # 
sim_opt         : smallint unsigned      # 
---
decoder="nnclassRaw"        : enum('nnclassRaw')            # 
ntrials=20                  : smallint                      # 
normalize=1                 : smallint                      # 
th_nonlin="rectify"         : enum('rectify','absolute')    # 
process="yes"               : enum('no','yes')              # m) do it or not
%}


classdef StatsSimParams < dj.Relvar
	methods

		function self = StatsSimParams(varargin)
			self.restrict(varargin{:})
		end
	end

end