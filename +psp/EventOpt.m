%{
# Optional parameters for event detection
psp_opt                      : smallint unsigned             # event detection option index
---
low_pass = 100               : float  # high pass filter
high_pass = 5                : float  # low pass filter
min_amp = 4                  : float  # minimum allowable amplitude for alpha functions in std of samples
max_amp= 0.0000000002        : float  # maximum allowable amplitude
min_tau= 5                   : float  # minimum allowable tau for alpha functions (in units of samples)
max_tau= 1000                : float  # maximum allowable tau
min_y_offset= -0.0000000001  : float  # minimum allowable yOffset for alpha functions (in units of mV)
max_y_offset= 0.0000000002   : float  # maximum allowable yOffset
min_decay= 5                 : float  # minimum allowable decay tau
max_decay= 1000              : float  # maximum allowable decay tau            
der_thresh= 0.000000000004   : float  # threshold used to determine if the change of derivative is great enough to separate out two alpha functions
closestepsps= 1000           : float  # second EPSP is not recognized if it is closer tha8 this value to the previous one (in units of samples)
errthresh= 0.004             : float  # threshold for standard error above which a fit with multiple alphas is attempted
datafiltertype= 3            : float  # 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay 
derfiltertype= 3             : float  # 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
datafilterlength= 100        : float  # length of data filter
derfilterlength= 100         : float  #  length of derivative filter
%}       
            
            
classdef EventOpt < dj.Lookup
end
