%{
# Patching experiments
session                      : timestamp           # recording timestamp
section                      : int                 # recoding subdivision
---
fs                           : float               # sampling rate
trial_duration               : float               # trial duration in seconds
trials                       : int                 # trial number
probes                       : int                 # probe number 
path                         : varchar(1024)       # file path
filename                     : varchar(255)        # file name
import_ts=CURRENT_TIMESTAMP  : timestamp           # don't edit
notes                        : varchar(1024)       # notes
%}       
            
            
classdef Exp < dj.Manual
    methods 
        function importData(self,filename)

            % get patch data
            [tree, data] = psp.ImportHEKAtoMat(filename);
        
            % fix data
            tree_map = ~cellfun(@isempty,tree);
            sub_idx = find(tree_map(:,3));
            tr = cell2mat(tree(tree_map(:,5),5));
            nprobes = length(unique([tr.TrSelfChannel]));

            % loop through sub-sections
            for isection = 1:length(sub_idx)
                
                % use correct data
                Data = data{1}((1:nprobes)+nprobes*(isection-1));
                
                % insert experiment key
                key = [];
                key.session = datestr((tree{1}.RoStartTimeMATLAB),'YYYY-mm-dd HH:MM:SS');
                key.section = isection;
                par_key = key;
                [key.path, key.filename] = fileparts(filename);
                key.fs = 1/tree{sub_idx(isection)+2,5}.TrXInterval;
                key.trial_duration = size(Data{1},1)/key.fs;
                key.trials = size(Data{1},2);
                key.probes = nprobes;
                fprintf('Importing section %d \n', isection)
                insert(psp.Exp,key)

                for iprobe = 1:nprobes
                    prob_key = par_key;
                    prob_key.probe = tree{sub_idx(isection)+iprobe+1,5}.TrSelfChannel;
                    prob_key.std = std(Data{iprobe}(:));
                    prob_key.mean = mean(Data{iprobe}(:));
                    prob_key.pipette_res = tree{sub_idx(isection)+iprobe+1,5}.TrPipetteResistance;
                    insert(psp.Probes,prob_key);

                     % loop through all traces
                    for itrial = 1:key.trials
                        tuple = par_key;
                        prob_key.probe = tree{sub_idx(isection)+iprobe+1,5}.TrSelfChannel;
                        tuple.trace = single(Data{iprobe}(:,itrial));
                        tuple.trial = tree{sub_idx(isection)+(itrial-1)*(nprobes+1)+1,4}.SwSweepCount;
                        tuple.probe = iprobe;
                        tuple.std = std(tuple.trace);
                        tuple.mean = mean(tuple.trace);
                        tuple.seal_res = tree{sub_idx(isection)+iprobe+(itrial-1)*(nprobes+1)+1,5}.TrSealResistance;
                        insert(psp.Traces,tuple);
                    end
                end
            end
        end
    end
end

