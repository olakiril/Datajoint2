%{
vis2p.MiceTrain (manual) # 
-> vis2p.Mice
start_timestamp : timestamp              # m) start timestamp of the trainning sesion
---
end_timestamp=null          : timestamp                     # m) end timestamp of the trainning sesion
train_type="ball"           : enum('ball,visual stimulus','ball,juice','handling','ball') # m) training type
train_notes                 : varchar(1023)                 # m) training notes
%}


classdef MiceTrain < dj.Relvar
	methods

		function self = MiceTrain(varargin)
			self.restrict(varargin{:})
		end
	end

end