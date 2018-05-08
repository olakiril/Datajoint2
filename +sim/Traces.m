%{
->sim.TraceParams       : smallint unsigned      #
->sim.RFParams
->stimulus.Movie
trace_id                          : smallint
---
trace                             : mediumblob                            #
%}


classdef Traces < dj.Computed
    properties
         keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.DecodeOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"'))))s
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
        end
    end
end

