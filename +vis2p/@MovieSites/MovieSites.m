%{
vis2p.MovieSites (computed) # 
-> vis2p.TracesGroup
-> vis2p.MovieParams
pThr            : float(6,0)             # 
snrThr          : smallint               # 
binsize         : smallint unsigned      # ms binsize for correlation computations
---
cor                         : mediumblob                    # c) correlation at  bin size
actual_binsize              : float                         # c) (ms) actual binsize used
synchrony                   : mediumblob                    # c) neural synchrony
pop_sparse                  : mediumblob                    # c) population sparseness
stim_length                 : double                        # c) the length of each movie in relative trace
act_sparse                  : mediumblob                    # c) tolhurst activity sparseness
kurtosis                    : mediumblob                    # c) kurtosis of the distribution
brain_state_covars_bs=CURRENT_TIMESTAMP: timestamp          # c) automatic timestamp, do not compute
%}


classdef MovieSites < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.TracesGroup*vis2p.MovieParams
	end

	methods(Access=protected)

		makeTuples( obj, key )
    end
    
    methods
		function self = MovieSites(varargin)
			self.restrict(varargin{:})
		end


	end

end