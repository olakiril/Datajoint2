%{
vis2p.ExtraTracesJake2 (computed) # 
-> vis2p.Scans
---
speed_trace=null            : mediumblob                    # ball speed
eye_trace=null              : mediumblob                    # eye movements [x y] (pixels)
pupil_trace=null            : mediumblob                    # pupil radius (pixels)
whisker_trace=null          : mediumblob                    # whisker movements [x y] (pixels)
shock_trace=null            : mediumblob                    # 
quality=null                : smallint                      # 
%}


classdef ExtraTracesJake2 < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Scans('problem_type = "none!" and state = "awake"').*vis2p.Experiments('process = "yes"').*vis2p.Movies.*vis2p.VisStims
	end

	methods(Access=protected)

        		makeTuples(obj,key)
    end
    methods
		function self = ExtraTracesJake2(varargin)
			self.restrict(varargin{:})
		end

		eye_dir = findEye(obj,key,eye_path)



		times = syncAffine(obj,tpObj,H5File,flipfreq)

		timeCorrector = syncVectors(obj,tpObj,H5File)

	end

end