%{
vis2p.LiquidDelivery (manual) # 
-> vis2p.BehSession
timestamp : bigint                 # timestamp of event
---
pulse_time      : smallint               # duration of pulse
microl_per_pulse: float                  # microliters per pulse time
liquid_type     : enum                   #water/juice
%}

classdef LiquidDelivery < dj.Relvar
     methods
        function self = LiquidDelivery(varargin)
           self.restrict(varargin{:})
        end
    end
end