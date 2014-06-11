%{
vis2p.OriTracesParams (lookup) # 
ori_opt         : tinyint unsigned       # m) opt ori index
---
shuffle=1000                : mediumint                     # m) number of shuffles
resp_delay=200              : smallint                      # m) the dalay between stimulus and response (msec
resp_period                 : smallint                      # m) time of response pediod (ms)
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef OriTracesParams < dj.Relvar
	methods

		function self = OriTracesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end