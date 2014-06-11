%{
vis2p.StatAreaData (computed) # 
-> vis2p.StatsSitesParams
mouse_id        : smallint unsigned      # 
exp_date        : date                   # 
scan_idx        : smallint unsigned      # 
trace_opt       : smallint unsigned      # 
movie_type      : varchar(10)            # the type of movie shown
stim_idx        : tinyint                # 
area            : varchar(5)             # 
area_opt        : smallint               # 
---
time                        : float(5,2)                    # 
trials                      : smallint                      # 
neurons                     : smallint                      # 
sigcorr                     : float(10,5)                   # 
mean                        : float(10,5)                   # 
variance                    : float(10,5)                   # c) variance
pspars                      : float(10,5)                   # c) population sparseness
lspars                      : float(10,5)                   # c) lifetime sparseness
pMulti                      : double                        # c) mutual information
pDist                       : double                        # c) classification performance
reliability                 : double                        # 
%}


classdef StatAreaData < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.StatsSites*vis2p.StatArea
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatAreaData(varargin)
			self.restrict(varargin{:})
		end


	end

end