%{
vis2p.ImageBar (imported) # 
-> vis2p.Scans
---
amp                         : longblob                      # amplitude of the fft phase spectrum
ang                         : longblob                      # angle of the fft phase spectrum
vessels=null                : mediumblob                    # 
%}


classdef ImageBar < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Scans('aim = "Intrinsic Imaging" and scan_prog = "MpScan" and stim_engine = "other"')
	end

	methods(Access=protected)
makeTuples( obj, key )
    end
    methods
		function self = ImageBar(varargin)
			self.restrict(varargin{:})
		end

		

		 plot(obj,varargin)

	end

end