%{
beh.StimPeriods (manual) # 
-> beh.Session
timestamp       : bigint                 # 
---
period_type=null            : char(255)                     # 
%}


classdef StimPeriods < dj.Relvar
end