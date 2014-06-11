%{
vis2p.StatAreaDecodeParams (lookup) # 
dec_opt         : tinyint                # 
---
ncells=0                    : tinyint                       # 
nclasses=2                  : tinyint                       # 
type="normal"               : enum('var','cov','rand','normal') # 
repetitions=30              : mediumint                     # 
classifier="nnclassRaw"     : enum('SVM','fisher','nnclassRaw') # m) the classifier used
process="yes"               : enum('no','yes')              # m) compute or not compute
%}


classdef StatAreaDecodeParams < dj.Relvar
	methods

		function self = StatAreaDecodeParams(varargin)
			self.restrict(varargin{:})
		end
	end

end