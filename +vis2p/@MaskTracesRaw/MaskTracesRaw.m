%{
vis2p.MaskTracesRaw (imported) # 
-> vis2p.MaskGroupRaw
masknum         : mediumint              # 
---
mask_type                   : enum('site','cells','neuropil','astrocyte','ephys','red','SST','PV','neuron') # c) where the data are coming from
calcium_trace               : mediumblob                    # i) unfiltered flourescence traces
annulus_trace=null          : mediumblob                    # i) calcium trace in an annulus around the cell
ephys_trace=null            : mediumblob                    # electrophysiology trace
red_trace=null              : mediumblob                    # i) red channel trace from the cell soma
%}


classdef MaskTracesRaw < dj.Relvar
        
	methods

	    plotSite(obj,tracetype)
		
        makeTuples( obj, key,traces,traces2,coordinates,masknums)

		function self = MaskTracesRaw(varargin)
			self.restrict(varargin{:})
		end
	end

end