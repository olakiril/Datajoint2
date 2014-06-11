%{
vis2p.StatsFrame (computed) # 
-> vis2p.TracesGroup
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
movie_type      : varchar(10)            # the type of movie shown
binsize         : smallint unsigned      # m) ms binsize for correlation computations
delay           : smallint               # m) calcium response delay
---
actual_binsize              : float                         # c) (ms) actual binsize used
traces                      : mediumblob                    # c) traces [bin trial movie cell]
sim_traces                  : mediumblob                    # c) simulated traces [bin trial movie rf]
ori_traces                  : mediumblob                    # c) simulated traces ori [bin trial movie ori]
im_k                        : mediumblob                    # c) pixel kurtosis [bin trial movie]
im_m                        : mediumblob                    # c) pixel mean [bin trial movie]
im_s                        : mediumblob                    # c) pixel std [bin trial movie]
im_c                        : mediumblob                    # c) pixel pwz corr [bin trial movie]
im_b                        : mediumblob                    # c) pixel bispectrum [bin trial movie]
INDEX(binsize,delay)
%}


classdef StatsFrame < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.TracesGroup*vis2p.StatsParams*vis2p.StatsFrameParams('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key )
    end
    
    methods
		function self = StatsFrame(varargin)
			self.restrict(varargin{:})
		end

		data = getTraceStruct(obj, key )

		traces = getTraces(obj, key )


	end

end