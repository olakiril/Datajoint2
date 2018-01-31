%{
# Lick performance with GonoGo tasks
->beh.Session
day                   : int32                         # day starts at 12AM
---
p                     : longblob                      # performance
p_shuffle             : longblob                      # chance performance
train_groups          : longblob                      # train group info
test_groups           : longblob                      # test group info
trial_info            : longblob                      # trial info [index, clip #, movie, bin #]
%}

classdef LickPerformance < dj.Computed
    properties
             keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.DecodeOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"'))))
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
        end
    end
end