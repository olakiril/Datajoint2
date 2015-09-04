%{
vis2p.BehSession (manual) # 
-> vis2p.Mice
session_timestamp : timestamp                 # timestamp of session
---
trial_interval     : smallint                   # time between trials (min)
response_interval     : smallint                   # buffer time for response (ms)
response_period     : smallint                   # response period time (min)
bad_delay     : smallint                   # punishment delay (min)
exp_type     : enum('Images','Orientation','BW','Freerun')                   # experiment type
stimuli     : char                   # stimuli to be presented, comma separated
rewarder_stimuli     : char                   # rewarded stimuli comma separated
setup     : tinyint                   # probe #
%}

classdef BehSession < dj.Relvar
     methods
        function self = BehSession(varargin)
           self.restrict(varargin{:})
        end
    end
end