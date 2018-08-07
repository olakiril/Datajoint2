%{
# Movie repeats evaluation
-> obj.Repeats
-> fuse.ScanSetUnit
-> stimulus.Clip
---
r                    : float                      # reliability
r_shuffle            : float                      # chance
p_shuffle            : float                      # p value
%}

classdef RepeatsUnit < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			% self.insert(key)
		end
	end

end