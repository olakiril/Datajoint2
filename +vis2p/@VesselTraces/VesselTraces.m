%{
vis2p.VesselTraces (imported) # 
-> vis2p.Movies
vmasknum        : mediumint              # 
---
contrast                    : decimal(6,4)                  # 
sharpness                   : decimal(6,4)                  # 
fit_qual                    : decimal(6,4)                  # 
x_trace                     : mediumblob                    # 
y_trace                     : mediumblob                    # 
radii_trace                 : mediumblob                    # 
timestamps                  : mediumblob                    # 
%}


classdef VesselTraces < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Movies.*vis2p.Experiments('process = "yes"').*vis2p.VisStims('exp_type = "MoviesExperiment" or exp_type = "MouseStatExperiment"').*vis2p.Scans('lens>12 and scan_prog = "MpScan"')
	end

	methods(Access=protected)

		makeTuples( obj, key)
        
    end
    
    methods
		function self = VesselTraces(varargin)
			self.restrict(varargin{:})
		end

		[T, binsize,Tr,uniMovies] = getTraces(obj,varargin)


	end

end