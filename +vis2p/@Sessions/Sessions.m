%{
vis2p.Sessions (manual) # list of sessions$
-> vis2p.Subjects
setup           : tinyint unsigned       # setup number
session_start_time: bigint               # start session timestamp
---
session_stop_time=null      : bigint                        # end of session timestamp
experimenter                : enum('James','Alex','Mani','Allison','Tori','Jacob','Dimitri','Cathryn','Manolis','Dennis','George','Shan') # name of person running exp
session_path                : varchar(255)                  # path to the data
session_datetime=null       : datetime                      # readable format of session start
hammer=0                    : tinyint                       # 
recording_software="Acquisition2.0": enum('Acquisition2.0','Hammer','Blackrock','Neuralynx') # software used to record the data
%}


classdef Sessions < dj.Relvar
	methods

		function self = Sessions(varargin)
			self.restrict(varargin{:})
		end
	end

end