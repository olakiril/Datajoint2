%{
vis2p.OriParams (computed) # my newest table
# add primary key here
-----
# add additional attributes
%}

classdef OriParams < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.VisStims('exp_type = "GratingExperiment"')
    end
    

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = OriParams(varargin)
			self.restrict(varargin{:})
		end


	end

end