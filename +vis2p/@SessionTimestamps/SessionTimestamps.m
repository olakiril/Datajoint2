%{
vis2p.SessionTimestamps (manual) # events recorded by the timestamper$
-> vis2p.Sessions
channel         : tinyint unsigned       # channel that received a message
timestamper_time: bigint                 # real time on computer
---
count                       : int unsigned                  # hardware count
%}


classdef SessionTimestamps < dj.Relvar
	methods

		function self = SessionTimestamps(varargin)
			self.restrict(varargin{:})
		end
	end

end