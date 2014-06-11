%{
vis2p.StatsImages (computed) # 
-> vis2p.StatsImagesParams
mouse_id        : smallint unsigned      # 
exp_date        : date                   # 
scan_idx        : smallint unsigned      # 
movie_num       : tinyint unsigned       # the number of the movie shown
movie_type      : varchar(10)            # the type of movie shown
trial_trigger   : tinyint                # 
---
im_kurtosis=null            : mediumblob                    # c) kurtosis of the pixel distribution
im_mean=null                : mediumblob                    # c) mean of the pixel intensity distribution
im_std=null                 : mediumblob                    # c) std of the pixel intesity distribution
im_pwz_corr=null            : mediumblob                    # c) second order statistics
im_bispectrum=null          : mediumblob                    # c) third order statistics
im_pwz_dist=null            : mediumblob                    # c) second order statistics
im_pdist=null               : mediumblob                    # 
im_diff=null                : mediumblob                    # 
im_motion=null              : mediumblob                    # 
%}


classdef StatsImages < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.TracesGroup('trace_opt = 17')*vis2p.StatsParams*vis2p.StatsImagesParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
    end
    
    methods
		function self = StatsImages(varargin)
			self.restrict(varargin{:})
		end

		out = getImStats(obj,key,functions,varargin)


	end

end