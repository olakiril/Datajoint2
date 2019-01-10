%{
# 
-> `pipeline_meso`.`__scan_set__unit`
%}


classdef MatchedUnits < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			 self.insert(key)
		end
	end

end