%{
vis2p.RFFit (computed) # 
-> vis2p.RFMap
---
snr                         : float                         # c) snr between rf and periphery 
gauss_fit                   : mediumblob                    # c) gauss fitting parameters
%}


classdef RFFit < dj.Relvar & dj.AutoPopulate

	properties(Constant)
		popRel = vis2p.RFMap('rf_opt_num =3')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = RFFit(varargin)
			self.restrict(varargin{:})
		end


		plot(obj,varargin)

	end

end