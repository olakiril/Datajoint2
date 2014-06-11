%{
vis2p.VisStims (imported) # 
-> vis2p.Scans
stim_idx        : tinyint unsigned       # the index of the stim file coupled to the scan file
---
stim_filename               : varchar(256)                  # i) the stimulus filename
stim_file                   : mediumblob                    # i) stim file without the swaptimes
exp_type="MultDimExperiment": enum('RodentGrating','MatrixMappingExperiment','MouseStatExperiment','DotMappingExperiment','MouseMovieExperiment','GratingExperiment','MoviesExperiment','BarMappingExperiment','other','MouseDotMapping','CenterSurround','MultDimExperiment') # i) 'RodentGrating' or 'MultDimExperiment'
frame_timestamps            : mediumblob                    # i) (ms) frame times synchronized to Labview clock
process="yes"               : enum('no','yes')              # m) process or not
INDEX(stim_idx)
%}


classdef VisStims < dj.Relvar & dj.AutoPopulate

	properties
		popRel = vis2p.Scans('stim_engine = "VisStim"  or stim_engine = "State" and problem_type = "none!" or stim_engine = "other" and problem_type = "none!"').*Experiments('process = "yes"')
	end

	methods(Access=protected)

		makeTuples( obj, key)
        
    end
    
    methods
		function self = VisStims(varargin)
			self.restrict(varargin{:})
		end

		gui( visStims )

		plot(obj)

	end

end