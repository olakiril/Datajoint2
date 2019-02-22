%{
# Patching experiments
parent                       : int                   # filename
---
fs                           : float               # sampling rate
trial_duration               : float               # trial duration in seconds
trials                       : int                 # trial number
probes                       : int                 # probe number 
path                         : int                 # file path
timestamp                    : timestamp           # recording timestamp
import_ts=CURRENT_TIMESTAMP  : timestamp           # don't edit
%}       
            
            
classdef Exp < dj.Manual
    methods 
        function importData(self,filename)

            [tree, data] = psp.ImportHEKAtoMat(filename);
            data = data{1};
            key = [];
            [key.path, key.parent] = fileparts(filename);
            
            % insert experiment key
            key.fs = 1/tree{5,5}.TrXInterval;
            key.trial_duration = size(data{1},1)/key.fs;
            key.trials = size(data{1},2);
            key.probes = size(data,2);
            key.timestamp = datestr((tree{1}.RoStartTimeMATLAB),'YYYY-mm-dd HH:MM:SS');
            insert(psp.Exp,key)
            
            for iprobe = 1:key.probes
                tuple = [];
                tuple.parent = key.parent;
                tuple.probe = tree{4+iprobe,5}.TrSelfChannel;
                tuple.std = std(data{iprobe}(:));
                tuple.mean = mean(data{iprobe}(:));
                tuple.pipette_res = tree{4+iprobe,5}.TrPipetteResistance;
                insert(psp.Probes,tuple);

                 % loop through all traces
                for itrial = 1:key.trials
                    tuple = [];
                    tuple.parent = key.parent;
                    tuple.trace = single(data{iprobe}(:,itrial));
                    tuple.trial = tree{4+(itrial-1)*(key.probes+1),4}.SwSweepCount;
                    tuple.probe = iprobe;
                    tuple.std = std(tuple.trace);
                    tuple.mean = mean(tuple.trace);
                    tuple.seal_res = tree{4+iprobe+(itrial-1)*(key.probes+1),5}.TrSealResistance;
                    insert(psp.Traces,tuple);
                end
            end
        end
    end
end

