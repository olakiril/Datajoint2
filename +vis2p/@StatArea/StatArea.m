%{
vis2p.StatArea (computed) #
-> vis2p.StatAreaParams
mouse_id        : smallint unsigned      #
exp_date        : date                   #
scan_idx        : smallint unsigned      #
---
v1o                         : mediumblob                    #
v1i                         : mediumblob                    #
v2o                         : mediumblob                    #
v2i                         : mediumblob                    #
fhMap                       : mediumblob                    #
reversal                    : mediumblob                    #
frev                        : mediumblob                    #
INDEX(mouse_id,exp_date,scan_idx)
INDEX(mouse_id,exp_date,scan_idx)
%}


classdef StatArea < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel = vis2p.MaskGroup*vis2p.StatAreaParams('process = "yes"').*vis2p.RFMapRaw  % !!! update the populate relation
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            %!!! compute missing fields for key here
            self.insert(key)
        end
        
    end
    
    methods
        function self = StatArea(varargin)
            self.restrict(varargin{:})
        end
    end
    
end