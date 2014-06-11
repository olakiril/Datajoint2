%{
vis2p.MovieTraces (computed) # 
-> vis2p.Traces
-> vis2p.MovieParams
binsize         : smallint unsigned      # ms binsize for correlation computations
undersample     : smallint               # 
---
outCorr                     : float                         # c)Correlations between different movietypes
inCorr                      : float                         # c) correlations between the same movie repetitions
outCorrP                    : float                         # c) significance of out correlation
inCorrP                     : float                         # c) significance of out correlation
mean_fr                     : float                         # c) the mean firing rate of the cell across one condition
ff                          : float                         # c) fano factor
cv                          : float                         # c) coefficient of variation
life_sparse                 : float                         # c) life sparseness
auto_corr                   : mediumblob                    # c)autocorrelation function (corr,lag,confidence)
kurtosis                    : float                         # c) kurtosis of the distribution
%}


classdef MovieTraces < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Traces*vis2p.MovieParams
    end
    
    

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = MovieTraces(varargin)
			self.restrict(varargin{:})
		end


	end

end