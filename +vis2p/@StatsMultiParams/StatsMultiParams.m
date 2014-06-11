%{
vis2p.StatsMultiParams (lookup) # 
multi_opt       : tinyint                # 
---
classifier="nnclassRaw"     : enum('nnclassRawSV','nnclassRaw') # m) the classifier used
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsMultiParams < dj.Relvar
	methods

		function self = StatsMultiParams(varargin)
			self.restrict(varargin{:})
		end
	end

end