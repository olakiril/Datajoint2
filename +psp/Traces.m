%{
# Patching traces
->psp.Exp
->psp.Probes
trial                       : int             # trial number
---
trace                       : mediumblob      # trace
std                         : float           # standard deviation of trace
%}      
            
classdef Traces < dj.Manual
end

