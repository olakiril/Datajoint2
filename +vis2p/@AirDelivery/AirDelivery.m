%{
vis2p.AirDelivery (manual) # 
-> vis2p.BehSession
timestamp : bigint                 # timestamp of event
---
pulse_time      : smallint               # duration of pulse
%}

classdef AirDelivery < dj.Relvar
     methods
        function self = AirDelivery(varargin)
           self.restrict(varargin{:})
        end
    end
end