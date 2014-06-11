%{
vis2p.TberPulses (manual) # trial sync pulses$
-> vis2p.Stimulation
tber_pulse_time : bigint                 # pulse time
---
%}


classdef TberPulses < dj.Relvar
	methods

		function self = TberPulses(varargin)
			self.restrict(varargin{:})
		end
	end

end