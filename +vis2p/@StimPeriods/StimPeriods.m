%{
vis2p.StimPeriods (manual) # 
-> vis2p.BehSession
timestamp : bigint                 # start timestamp of period
---
period_type     : char                   # discription of period
%}

classdef StimPeriods < dj.Relvar
     methods
        function self = StimPeriods(varargin)
           self.restrict(varargin{:})
        end
    end
end