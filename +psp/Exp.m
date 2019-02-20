%{
# Patching experiments
parent                       : varchar             # filename
---
fs                           : float               # sampling rate
ts=CURRENT_TIMESTAMP         : timestamp                     # don't edit
%}       
            
            
classdef Exp < dj.Lookup
    methods
        function insert(filename)
            l = load(filename);
            names = fieldnames(l);
            trials = cellfun(@(x) str2double(x{1}{1}),regexp(names, 'Trace_1_1_(\d*)_\d*','tokens'));
            probes = cellfun(@(x) str2double(x{1}{1}),regexp(names, 'Trace_1_1_\d*_(\d*)','tokens'));
            for itrace = 1:length(l)
            end
        end
    end
end

