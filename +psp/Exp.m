%{
# Patching experiments
parent                       : int                   # filename
---
fs                           : float               # sampling rate
trial_duration               : float               # trial duration in seconds
trials                       : int                 # trial number
probes                       : int                 # probe number 
path                         : int            # file path
ts=CURRENT_TIMESTAMP         : timestamp           # don't edit
%}       
            
            
classdef Exp < dj.Manual
    methods 
        function importData(self,filename)
            data = load(filename);
            key = [];
            [key.path, key.parent] = fileparts(filename);
            names = fieldnames(data);
            trials = cellfun(@(x) str2double(x{1}{1}),regexp(names, 'Trace_1_1_(\d*)_\d*','tokens'));
            probes = cellfun(@(x) str2double(x{1}{1}),regexp(names, 'Trace_1_1_\d*_(\d*)','tokens'));
            data = struct2cell(data);
            un_probes = unique(probes);
            
            % insert experiment key
            key.fs = 1/median(diff(data{1}(:,1)));
            key.trial_duration = length(data{1}(:,2))/key.fs;
            key.trials = length(unique(trials));
            key.probes = length(un_probes);
            insert(psp.Exp,key)
            
            for probe = un_probes'
                tuple = [];
                tuple.parent = key.parent;
                tuple.probe = probe;
                tuple.std = std(cell2mat(cellfun(@(x) x(:,2),data(un_probes==probe),'uni',0)));
                insert(psp.Probes,tuple);
            end
            
            % loop through all traces
            for itrace = 1:length(trials)
                tuple = [];
                tuple.parent = key.parent;
                tuple.trace = single(data{itrace}(:,2));
                tuple.trial = trials(itrace);
                tuple.probe = probes(itrace);
                tuple.std = std(tuple.trace);
                insert(psp.Traces,tuple);
            end
        end
    end
end

