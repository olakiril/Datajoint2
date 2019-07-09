%{
# Movie repeats evaluation
-> simf.Repeats
-> simf.ActivityTrace
---
r                    : float                      # reliability
%}

classdef RepeatsUnit < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			% self.insert(key)
		end
	end

end