%{
mov3d.ReduceDMOpt (lookup) # 
red_opt       : smallint unsigned      # 
---
brief="fill out"            : varchar(127)                        # short description, to be displayed in menus
binsize=500                 : float                               # time window in ms to compute the response
reduce_method="tsne"        : enum('tsne')                        # dimensionality reduction method
dimensions=3                : mediumint                           # final dimensions
initial_dims=50             : mediumint                           # initial_dims after PCA
perplexity=30               : mediumint                           # perplexity
gauss_win=1500              : int                                 # trial time bluring in msec
process="yes"               : enum('no','yes')                    # do it or not
%}


classdef ReduceDMOpt < dj.Relvar
	methods
		function self = DecodeTimeOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end