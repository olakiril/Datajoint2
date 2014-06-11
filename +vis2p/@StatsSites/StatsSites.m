%{
vis2p.StatsSites (computed) # 
-> vis2p.TracesGroup
-> vis2p.StatsSitesParams
movie_type      : varchar(10)            # the type of movie shown
stim_idx        : tinyint                # 
trial_trigger   : tinyint                # 
---
time                        : float(5,2)                    # 
trials                      : smallint                      # 
neurons                     : smallint                      # 
synchrony                   : float(10,5)                   # c) neural synchrony
eudist                      : float(10,5)                   # c)average eucleadian distance
sigcorr                     : float(10,5)                   # 
mean                        : float(10,5)                   # 
variance                    : float(10,5)                   # c) variance
pspars                      : float(10,5)                   # c) population sparseness
pzero                       : float(10,5)                   # c) probability of zero neurons responding
pkurt                       : float(10,5)                   # c) kurtosis
lspars                      : float(10,5)                   # c) lifetime sparseness
lzero                       : float(10,5)                   # c) probability of zero response
lkurt                       : float(10,5)                   # c) lifetime kurtosis
pspars_tr                   : float(10,5)                   # c) population sparseness
pzero_tr                    : float(10,5)                   # c) probability of zero neurons responding
pkurt_tr                    : float(10,5)                   # c) kurtosis
lspars_tr                   : float(10,5)                   # c) lifetime sparseness
lzero_tr                    : float(10,5)                   # c) probability of zero response
lkurt_tr                    : float(10,5)                   # c) lifetime kurtosis
varexp                      : float(10,5)                   # 
INDEX(mouse_id,exp_date,scan_idx,movie_type)
INDEX(mouse_id,exp_date,scan_idx,stats_opt,trace_opt,movie_type)
INDEX(movie_type)
%}


classdef StatsSites < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.TracesGroup*vis2p.StatsParams*vis2p.StatsSitesParams('process = "yes"').*vis2p.Scans('problem_type = "none!"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
        
    end
    
    methods
		function self = StatsSites(varargin)
			self.restrict(varargin{:})
		end

		traces = fetchTraces(obj)

		[T, binsize] = getExtraTraces(obj,type,varargin)

		[T, binsize] = getExtraTracesJake(obj,type,varargin)

		[T, binsize] = getSpikes(obj,varargin)

		out = getStimulus(obj,key,varargin)

		out = getStimulusTrace( obj, key )

		[T, binsize,Tr,uniMovies] = getTraces(obj,varargin)

		imTrials(obj,varargin)

		plot( obj )

		plot2(obj)

		plot3(obj,varargin)

		plotSingle(obj,varargin)

		plotTrials(obj )

		plotnum(obj,varargin)

		trialCorr(obj )

	end

end