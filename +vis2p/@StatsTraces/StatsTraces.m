%{
vis2p.StatsTraces (computed) # 
-> vis2p.Traces
-> vis2p.StatsParams
-> vis2p.StatsTracesParams
---
outCorr                     : float                         # c)Correlations between different movietypes
inCorr                      : float                         # c) correlations between the same movie repetitions
outCorrP                    : float                         # c) significance of out correlation
inCorrP                     : float                         # c) significance of out correlation
mean                        : mediumblob                    # c) the mean firing rate of the cell across one condition
life_sparse                 : float                         # c) life sparseness
auto_corr                   : mediumblob                    # c)autocorrelation function (corr,lag,confidence)
kurtosis                    : float                         # c) kurtosis of the distribution
variance                    : mediumblob                    # c) variance of the trials
life_sparse_avg             : float                         # c) life sparseness
pzero                       : float                         # m) probability of zero response
pzero_avg                   : float                         # m) probability of zero response, on averaged trials
kurtosis_avg                : float                         # c) kurtosis of the avg distribution
INDEX(binsize,undersample)
INDEX(mouse_id,exp_date,scan_idx,stim_idx,movie_type)
%}


classdef StatsTraces < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Traces*vis2p.StatsParams*vis2p.StatsTracesParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsTraces(varargin)
			self.restrict(varargin{:})
		end


		plot(obj )

		plot2(obj,amp)

		plotSingle(obj,cells)

	end

end