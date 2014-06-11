%{
vis2p.MoviePresents (computed) # 
-> vis2p.MovieParams
repeat_num      : mediumint unsigned     # 
---
movie_start_time            : mediumint unsigned            # i) the number of the movie shown in this trial
movie_times                 : mediumblob                    # c) movieframe timestamps in win time 
%}


classdef MoviePresents < dj.Relvar
	methods

	makeTuples( obj, key, stim)


		function self = MoviePresents(varargin)
			self.restrict(varargin{:})
		end
	end

end