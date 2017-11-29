%{
mov3d.RFFilterOpt (lookup) # 
rf_opt       : smallint unsigned      # 
---
brief="fill out"            : varchar(127)                        # short description, to be displayed in menus
binsize=500                 : float                               # time window in ms to compute the response
rf_thr=2                    : float                               # standard deviation of population RF
process="yes"               : enum('no','yes')                    # do it or not
%}


classdef RFFilterOpt < dj.Relvar
	methods
		function self = RFFilterOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end