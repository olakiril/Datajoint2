%{
vis2p.MaskGroupRaw (imported) #
-> vis2p.Movies
---
frame_timestamps            : mediumblob                    # i) (ms) frame times synchronized to Labview clock
fps                         : double                        # i) frames per second
%}


classdef MaskGroupRaw < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel = vis2p.Movies.*vis2p.Scans('scan_prog = "AOD"').*vis2p.Experiments('process = "yes"').*vis2p.VisStims
    end
    
    methods(Access=protected)
        makeTuples( obj, key )
    end
    methods
        function self = MaskGroupRaw(varargin)
            self.restrict(varargin{:})
        end
    end
    
end