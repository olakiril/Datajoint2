%{
vis2p.StatsSimTracesParams (lookup) # 
sim_traces_opt  : smallint unsigned      # 
---
movie_resize=1.0000         : float(5,4)                    # m) image resize ratio
filter_type="sparsenet"     : enum('Pixels','DoG,Pixels','DoG,sparsenet','sparsenet') # 
frame_step=3                : smallint                      # 
movie_path="Q:/MouseMovie/MacNew/": varchar(512)            # 
filter_size=0.5219          : float(5,4)                    # m) filter relative size to the image
process="yes"               : enum('no','yes')              # m) do it or not
%}


classdef StatsSimTracesParams < dj.Relvar
	methods

		function self = StatsSimTracesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end