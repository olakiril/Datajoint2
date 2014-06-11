%{
vis2p.RFPresents (computed) # 
-> vis2p.RFParams
repeat_num      : mediumint unsigned     # 
---
dot_time                    : bigint                        # i) timestamp of the stimulus on labview times
dot_loc_x                   : smallint unsigned             # i) x location in screen
dot_loc_y                   : smallint                      # i) x location in screen
dot_color                   : smallint unsigned             # i) color of the stimulus
%}


classdef RFPresents < dj.Relvar
	methods

		makeTuples( obj, key, stim, indx )

		function self = RFPresents(varargin)
			self.restrict(varargin{:})
		end
	end

end