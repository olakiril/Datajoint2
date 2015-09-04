%{
vis2p.OptImageBar (imported) # 
-> vis2p.Scans
---
amp                         : longblob                      # amplitude of the fft phase spectrum
ang                         : longblob                      # angle of the fft phase spectrum
vessels=null                : mediumblob                    # 
direction=0                 : smallint                      # 
%}


classdef OptImageBar < dj.Relvar & dj.AutoPopulate

	properties (Constant)
		popRel = vis2p.Scans('problem_type = "none!"')...
            *vis2p.VisStims('exp_type = "BarMappingExperiment" or exp_type = "FancyBarExperiment"')
    end

    
    methods(Access=protected)
        
		makeTuples( obj, key )
        
    end
	methods

		function self = OptImageBar(varargin)
			self.restrict(varargin{:})
		end


		 [iH, iS, iV] = plot(obj,varargin)

	end

end