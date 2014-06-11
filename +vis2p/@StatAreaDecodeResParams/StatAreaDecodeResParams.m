%{
vis2p.StatAreaDecodeResParams (lookup) # 
dec_opt_res     : tinyint                # 
---
area="V1"                   : enum('V1','V2')               # 
area_classifier="nnclassRaw": enum('SVM','fisher','nnclassRaw') # 
trials="all"                : enum('true','false','all')    # 
source="Stim"               : enum('Stim','V2','V1')        # 
source_type="normal"        : enum('cov','rand','normal')   # 
process="yes"               : enum('no','yes')              # m) compute or not compute
discription=null            : varchar(255)                  # 
%}


classdef StatAreaDecodeResParams < dj.Relvar
	methods

		function self = StatAreaDecodeResParams(varargin)
			self.restrict(varargin{:})
		end
	end

end