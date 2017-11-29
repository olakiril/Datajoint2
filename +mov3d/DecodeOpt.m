%{
mov3d.DecodeOpt (lookup) # 
dec_opt       : smallint unsigned      # 
---
brief="fill out"            : varchar(127)                        # short description, to be displayed in menus
binsize=500                 : float                               # time window in ms to compute the response
decode_method="nnclassRaw"  : enum('nnclassRaw','nnclassRawSV','nnclass')   # decoding method
select_method="all"         : enum('all','subsample','expand')    # cell selection 
trial_bins=1                : mediumint                           # trial grouping
trial_method="random"       : enum('random','sequential')         # trial selection method
process="yes"               : enum('no','yes')                    # do it or not
rf_opt=0                    : smallint unsigned                   # restrict anaysis to population RF
chance=0                    : tinyint                             # chance performance
%}


classdef DecodeOpt < dj.Relvar
	methods
		function self = DecodeOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end