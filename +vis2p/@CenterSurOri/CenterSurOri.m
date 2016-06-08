%{
vis2p.CenterSurOri (computed) #
-> vis2p.Traces
-> vis2p.CenterSurParams
---
oriIn                       : mediumblob                    # c) the raw data [trials,oris]
PotiIn                      : float(4,3)                    # c) the significance of tuning from oti
PdmIn                       : float(4,3)                    # c) oreffered direction of motion
oriOut                      : mediumblob                    # c) the raw data [trials,oris]
PotiOut                     : float(4,3)                    # c) the significance of tuning from oti
PdmOut                      : float(4,3)                    # c) oreffered direction of motion
maxTraces                   : mediumblob                    #
INDEX(mouse_id,exp_date,scan_idx)
INDEX(mouse_id,exp_date,scan_idx,center_sur_opt,trace_opt)
%}


classdef CenterSurOri < dj.Relvar & dj.AutoPopulate
    
    properties(Constant)
        popRel = vis2p.Traces*vis2p.CenterSurParams('process = "yes"') & (vis2p.Scans('problem_type = "none!"') & vis2p.CenterSurTrials('repeat_num = 1'))
    end
    
    methods
        
        function self = CenterSurOri(varargin)
            self.restrict(varargin{:})
        end
        
        [T,oris, binsize] = getTraces(obj,varargin)
        
        plot(obj,type)
        
    end
    methods(Access=protected)
        makeTuples( obj, key ) 
    end
    
end