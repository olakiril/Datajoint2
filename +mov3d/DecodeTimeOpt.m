%{
mov3d.DecodeTimeOpt (lookup) # 
dec_opt       : smallint unsigned      # 
---
brief="fill out"            : varchar(127)                        # short description, to be displayed in menus
binsize=500                 : float                               # time window in ms to compute the response
decode_method="nnclassRaw"  : enum('nnclassRaw','nnclassRawSV','nnclass')   # decoding method
trial_bins=1                : mediumint                           # trial grouping
trial_method="random"       : enum('random','sequential')         # trial selection method
process="yes"               : enum('no','yes')                    # do it or not
%}


classdef DecodeTimeOpt < dj.Relvar
	methods
		function self = DecodeTimeOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end