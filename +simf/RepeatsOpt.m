%{
# Repeats options table
rep_opt                           : smallint unsigned      # 
---
brief="fill out"                  : varchar(127)                        # short description, to be displayed in menus
method="explainedVar"             : enum('explainedVar','corr','oracle')# reliability method
process="yes"                     : enum('no','yes')                    # do it or not
sample_sz=1000                    : int                            # bootstrap shuffling number
noise="@(x) wblrnd(x,0.5)/2"      :varchar(1024)                        # noise type
%}


classdef RepeatsOpt < dj.Lookup
end