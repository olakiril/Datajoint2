%{
vis2p.StatsPresents (computed) # 
-> vis2p.StatsParams
repeat_num      : mediumint unsigned     # 
---
movie_num                   : mediumint unsigned            # i) the number of the movie shown in this trial
movie_times                 : mediumblob                    # c) movieframe timestamps in win time 
%}


classdef StatsPresents < dj.Relvar
	methods

        makeTuples( obj, key, stim, movieStat,trialTrigger )


		function self = StatsPresents(varargin)
			self.restrict(varargin{:})
		end
	end

end