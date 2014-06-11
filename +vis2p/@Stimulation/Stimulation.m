%{
vis2p.Stimulation (manual) # visual stimulations$
subject_id      : int unsigned           # unique identifier for subject
setup           : tinyint unsigned       # setup number
session_start_time: bigint               # start session timestamp
stim_start_time : bigint                 # timestamp for stimulation start
---
stim_stop_time=null         : bigint                        # end of stimulation timestamp
stim_path                   : varchar(255)                  # path to the stimulation data
exp_type                    : varchar(255)                  # type of experiment
total_trials=null           : int unsigned                  # total number of trials completed
correct_trials=null         : int unsigned                  # number of correct trials
incorrect_trials=null       : int unsigned                  # number of incorrect trials
%}


classdef Stimulation < dj.Relvar
	methods

		function self = Stimulation(varargin)
			self.restrict(varargin{:})
		end
	end

end