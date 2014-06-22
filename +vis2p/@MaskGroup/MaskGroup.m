%{
vis2p.MaskGroup (imported) #
-> vis2p.Movies
---
%}


classdef MaskGroup < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel = vis2p.Movies.*vis2p.Scans('lens>12 or scan_prog = "Unirec"').*vis2p.Experiments('process = "yes"')
    end
    
    methods(Access=protected)
        makeTuples( obj, key )
    end
    methods
        function self = MaskGroup(varargin)
            self.restrict(varargin{:})
        end
        
        insertCell( obj, key,masktype)
    end
    
end