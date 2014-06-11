%{
vis2p.Stacks (imported) # 
-> vis2p.Scans
---
xsize                       : smallint unsigned             # i) the x frame size of the scan
ysize                       : smallint unsigned             # i) the y frame size of the scan
interval                    : smallint                      # the interval in z, positive means deeper direction
steps                       : smallint                      # number of steps in z
averages                    : smallint                      #  number of averages per xy plane
%}


classdef Stacks < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Scans('aim = "stack" and problem_type = "none!"').*vis2p.Experiments('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = Stacks(varargin)
			self.restrict(varargin{:})
		end


		tpr = tpReader( obj )

	end

end