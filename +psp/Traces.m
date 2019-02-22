%{
# Patching traces
->psp.Exp
->psp.Probes
trial                       : int             # trial number
---
trace                       : mediumblob      # trace
std                         : float           # standard deviation of trace
mean                        : float           # mean of trace
seal_res                    : float           # seal resistance
%}      
            
classdef Traces < dj.Manual
end

