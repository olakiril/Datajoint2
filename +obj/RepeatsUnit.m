%{
# Movie repeats evaluation
-> obj.Repeats
-> fuse.ScanSetUnit
-> stimulus.Clip
---
r                    : longblob                      # reliability
r_shuffle            : longblob                      # chance
p_shuffle            : longblob                      # p value
%}

classdef RepeatsUnit < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			% self.insert(key)
		end
	end

end