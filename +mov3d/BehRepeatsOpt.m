%{
mov3d.BehRepeatsOpt (lookup) # 
rep_opt       : smallint unsigned      # 
---
brief="fill out"                  : varchar(127)                        # short description, to be displayed in menus
binsize=500                       : float                               # time window in ms to compute the response
method="explainedVar"             : enum('explainedVar','corr')         # reliability method
process="yes"                     : enum('no','yes')                    # do it or not
restrict_rf=0                     : float                               # restrict anaysis to population RF
beh_thr=0                         : float                               # standard deviations from mean running threshold
%}


classdef BehRepeatsOpt < dj.Relvar
	methods
		function self = BehRepeatsOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end