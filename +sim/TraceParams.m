%{
trace_opt       : smallint unsigned      #
---
bin                               : smallint                            # temporal window (ms)
nonlinearity                      : enum('thr_exp')                     # nonlinearity type
brief="fill out"                  : varchar(127)                        # short description, to be displayed in menus
process="yes"                     : enum('no','yes')                    # do it or not
%}


classdef TraceParams < dj.Lookup
end

