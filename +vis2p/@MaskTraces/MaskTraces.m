%{
vis2p.MaskTraces (imported) # 
-> vis2p.MaskGroup
masknum         : mediumint              # 
---
mask_type                   : enum('site','cells','neuropil','astrocyte','ephys','red','SST','PV','neuron') # c) where the data are coming from
calcium_trace               : mediumblob                    # i) unfiltered flourescence traces
annulus_trace=null          : mediumblob                    # i) calcium trace in an annulus around the cell
ephys_trace=null            : mediumblob                    # electrophysiology trace
red_trace=null              : mediumblob                    # i) red channel trace from the cell soma
%}


classdef MaskTraces < dj.Relvar
	methods

	    plotSite(obj,tracetype)
		
          makeTuples( obj, key, cellFinder, type,shareMask)


		function self = MaskTraces(varargin)
			self.restrict(varargin{:})
		end
	end

end