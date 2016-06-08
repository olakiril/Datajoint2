%{
olf.StimPeriods (manual) # 
-> olf.Session
timestamp                      : bigint        # time from the start of the session in ms
---
trial=null                     : int           # trial number
stimulus=null                  : varchar       # Stimuli indexes that were used
%}


classdef StimPeriods < dj.Relvar
    methods
        function self = StimPeriods(varargin)
            self.restrict(varargin{:})
        end
    end
end