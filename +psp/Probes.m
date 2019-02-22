%{
# Patching probes
->psp.Exp
probe                       : int             # probe number
---
std                         : float           # standard deviation of probe
mean                         : float           # mean of probe
pipette_res                   : float           # pipette resistance
%}      
            
classdef Probes < dj.Manual
end

