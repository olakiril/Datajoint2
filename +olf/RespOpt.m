%{
olf.RespOpt (lookup) # 
resp_opt       : smallint unsigned      # 
---
brief="fill out"            : varchar(127)                  # short description, to be displayed in menus
baseline_period             : double                        # Period before onset to average baseline signal (ms)
response_period             : double                        # Period to average signal (ms)
process="yes"               : enum('no','yes')              # do it or not
%}


classdef RespOpt < dj.Relvar
	methods

		function self = RespOpt(varargin)
			self.restrict(varargin{:})
		end
	end

end