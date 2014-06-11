%{
vis2p.StatsAwakeParams (lookup) # 
awake_opt       : smallint               # 
---
trace_opt                   : smallint                      # 
stats_opt                   : smallint                      # 
ball_thr=null               : double(5,3) unsigned          # m) ball movement threshold (?)
whisker_thr=null            : double(5,1) unsigned          # m) whisker movement threshold (arb)
eye_thr=null                : double(5,2) unsigned          # m) eye movement threshold (pixels)
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsAwakeParams < dj.Relvar
	methods

		function self = StatsAwakeParams(varargin)
			self.restrict(varargin{:})
		end
	end

end