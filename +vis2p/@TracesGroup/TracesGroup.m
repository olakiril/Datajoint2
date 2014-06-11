%{
vis2p.TracesGroup (computed) # 
-> vis2p.MaskGroup
-> vis2p.TracesOpt
---
%}


classdef TracesGroup < dj.Relvar & dj.AutoPopulate

	properties
		popRel = (vis2p.MaskGroup*vis2p.TracesOpt('process = "yes"')).*vis2p.Experiments('process = "yes"')
	end

	methods(Access=protected)

		makeTuples(obj,key)
        
    end
    
    methods
		function self = TracesGroup(varargin)
			self.restrict(varargin{:})
		end


		plotSite(obj)

		plotTraces(obj,key,trace_opt,colors)

	end

end