%{
olf.Session (manual) #
-> olf.Mice
session_timestamp: timestamp             #
---
filename=null               : varchar                       # the name of the scan file
depth=null                  : int                           # depth of the recording in microns
pulse_duration=null         : smallint                      # odor presentation time (ms)
intertrial_interval=null    : smallint                      # time between odors (ms)
setup=null                  : tinyint                       # setup number
lens=null                   : tinyint                       # lens magnification
notes                       : varchar(256)                  # optional free-form comments
path                        : varchar(255)                  # path to the file
%}


classdef Session < dj.Relvar
    
    methods
        function self = Session(varargin)
            self.restrict(varargin{:})
        end
    end
end