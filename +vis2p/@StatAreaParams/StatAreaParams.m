%{
vis2p.StatAreaParams (lookup) # 
area_opt        : smallint               # param index
---
gwin=200                    : smallint unsigned             # m) gaussian window size to convolve the RF Site Map
glim=50                     : smallint                      # m) distance from reversal in pixels difining inner bounds
process="yes"               : enum('no','yes')              # m) compute or not compute
INDEX(gwin)
%}


classdef StatAreaParams < dj.Relvar
	methods

		function self = StatAreaParams(varargin)
			self.restrict(varargin{:})
		end
	end

end