%{
vis2p.Licks (manual) # 
-> vis2p.Mice
timestamp : bigint                 # timestamp of event
---
session_timestamp : timestamp                 # timestamp of session
exp_date        : date                   # experiment date
%}

classdef Licks < dj.Relvar
     methods
        function self = Licks(varargin)
           self.restrict(varargin{:})
        end
    end
end