%{
olf.MaskTraces (imported) # 
-> olf.Scan
masknum         : mediumint              # 
---
mask                        : blob                          # mask used to extract traces
mask_type                   : enum('neuropil','astrocyte','neuron') # c) where the data are coming from
trace_ch1                   : mediumblob                    # unfiltered flourescence traces from channel 1
trace_ch2                   : mediumblob                    # unfiltered flourescence traces from channel 2
xpos                        : int                           # cell center x position
ypos                        : int                           # cell center y position
zpos                        : int                           # cell center z position
radius                      : double                        # radius of the cell in microns
%}


classdef MaskTraces < dj.Relvar
	methods
		
        function makeTuples( obj, key, masknum, mask, mask_type, trace_ch1, trace_ch2,xpos,ypos,zpos,rad)
            tuple = key;
            tuple.mask = mask;
            tuple.masknum = masknum;
            tuple.mask_type = mask_type;
            tuple.trace_ch1 = trace_ch1;
            tuple.trace_ch2 = trace_ch2;
            tuple.xpos = xpos;
            tuple.ypos = ypos;
            tuple.zpos = zpos;
            tuple.radius = rad;
            insert( obj, tuple );
        end

		function self = MaskTraces(varargin)
			self.restrict(varargin{:})
		end
	end

end