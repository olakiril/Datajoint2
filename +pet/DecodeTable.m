%{
# 
-> `pipeline_stimulus`.`_sync`
-> `pipeline_fuse`.`__scan_done`
-> pet.Areas
-> pet.DecodeOpts
---
p                           : longblob                      # performance
p_shuffle                   : longblob                      # chance performance
train_groups                : longblob                      # train group info
test_groups                 : longblob                      # test group info
trial_info                  : longblob                      # trial info [index, clip #, movie, bin #]
score                       : longblob                      # svm score (distance to boundary)
%}


classdef DecodeTable < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			 self.insert(key)
		end
	end

end