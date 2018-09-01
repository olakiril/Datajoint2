%{
#
resp_opt             : int
---
high_pass_fr             : double   #
low_pass_fr              : double   #
peak_min_interval        : double   #
trough_min_interval      : double   #
peak_thr_prctile         : double   #
trough_thr_prctile       : double   #
peak_thr_factor          : double   #
trough_thr_factor        : double   #
peak_height_thr_prctile  : double   #
trough_height_thr_prctile: double   #
peak_height_thr_factor   : double   #
trough_height_thr_factor : double   #
%}


classdef RespiOpt < dj.Lookup
    methods (Access=protected)
    end
end