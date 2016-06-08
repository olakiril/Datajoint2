%{
olf.Stimuli (manual) # 
-> olf.Session
stimulus_index              : int                           # Unique stimulus identification number
---
stimulus_description        : varchar                       # stimulus description
%}

classdef Stimuli < dj.Relvar
     methods
        function self = Stimuli(varargin)
           self.restrict(varargin{:})
        end
    end
end