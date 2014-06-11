%{
vis2p.EphysTraces (imported) # 
-> vis2p.MaskGroup
---
masknum                     : smallint                      # m) mask number of recorded cell
spike_times                 : mediumblob                    # detected spikes
problem=null                : enum('none','bad signal','other cell') # 
%}


classdef EphysTraces < dj.Relvar & dj.AutoPopulate

	properties(Constant)
		popRel = (vis2p.MaskGroup.*vis2p.Scans('aim = "patching"')).*vis2p.Experiments('process = "yes"')
	end

	methods(Access=protected)

        		makeTuples( obj, key)
    end
    methods
		function self = EphysTraces(varargin)
			self.restrict(varargin{:})
		end

		[traces, keys] = getCaTraces(obj,varargin)



		out = plotCaTrace(obj,varargin)

		fhout = plotMask(obj,varargin)

		[regOUT, pOUT, rOUT] = regressCa(obj,varargin)

	end

end