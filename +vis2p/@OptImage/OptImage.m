%{
vis2p.OptImage (imported) # 
-> vis2p.Scans
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
---
spot_amp                    : longblob                      # (percent) response magnitudes for each spot
spot_r2                     : longblob                      # r-squared of total response
spot_fp                     : longblob                      # total response p-value (F-test)
spot_psth                   : longblob                      # average stimulus-locked response
vessels=null                : mediumblob                    # 
%}


classdef OptImage < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.VisStims('exp_type = "MultDimExperiment"').*vis2p.Scans('aim = "Intrinsic Imaging"')
	end

	methods(Access=protected)

        
		makeTuples( obj, key )
        
    end
    
    methods
		function self = OptImage(varargin)
			self.restrict(varargin{:})
		end


		 plot(obj,varargin)

	end

end