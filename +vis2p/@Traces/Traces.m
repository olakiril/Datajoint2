%{
vis2p.Traces (computed) #
-> vis2p.TracesGroup
masknum         : mediumint              #
---
trace                       : mediumblob                    # c) filtered traces using opt
ts=CURRENT_TIMESTAMP        : timestamp                     # c) extaction timestamp
quality=null                : float                         # m) correlation between reconstructed and raw trace
%}


classdef Traces < dj.Relvar & dj.AutoPopulate
    
    
    methods
        function self = Traces(varargin)
            self.restrict(varargin{:})
        end
        
        p = plot(obj,opts)
        
        p = plotStim(obj,opts,cind,lp,fig)
    end
    methods(Access=protected)
        
        makeTuples(self, key)
           
    end
    
end