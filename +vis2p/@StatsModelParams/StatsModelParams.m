%{
vis2p.StatsModelParams (lookup) # 
class_opt       : tinyint                # 
---
type="trace"                : enum('var','cov','rand','trace') #  type of trace computation
repetitions=10              : smallint                      # 
classes=999                 : smallint                      # m) the number of underclasses
classrep=1                  : smallint                      # m) the number of repetitions for underclassing
cells="all"                 : enum('vary','all')            # 
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatsModelParams < dj.Relvar
	methods

		function self = StatsModelParams(varargin)
			self.restrict(varargin{:})
		end
	end

end