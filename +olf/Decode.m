%{
olf.Decode (computed) # Odor decoding
-> olf.DecodeOpt
-> olf.RespOpt
-> preprocess.Spikes
-> preprocess.MaskType
---
p                     : longblob                      # performance
p_shuffle             : longblob                      # chance performance
train_groups          : longblob                      # train group info
test_groups           : longblob                      # test group info
%}

classdef Decode < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  =  (preprocess.Spikes * preprocess.MaskType  & preprocess.MaskClassificationMaskType) ...
            * pro(olf.RespOpt  & 'process = "yes"','process->pro') * (olf.DecodeOpt & 'process = "yes"') & olf.Sync
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            [decoder,k_fold,shuffle,train_set,test_set,repetitions] = ...
                fetch1(olf.DecodeOpt & key,...
                'decoder','k_fold','shuffle','train_set','test_set','repetitions');
            
            % get Data
            [resp_on, resp_off, stimuli] = fetchn(olf.OlfResponses * preprocess.MaskClassificationMaskType & key,'resp_on','resp_off','stimuli');
            resp_on = cell2mat(cellfun(@(x) permute(x,[3 1 2]),resp_on,'uni',0));
            resp_off = cell2mat(cellfun(@(x) permute(x,[3 1 2]),resp_off,'uni',0));
            resp = [resp_on;resp_off];
            Stims = stimuli{1};
            Traces = squeeze(num2cell(permute(resp,[1 3 2]),[1 2]))';

            [train_groups,test_groups] = getGroups(self,Stims,train_set,test_set);
            
            % check for correct stimuli
            cl = [train_groups{:}]; cl = [cl{:}];
            if ~isempty(test_groups)
                cl2 = [test_groups{:}]; cl2 = [cl2{:}];
            else
                cl2 = cl;
            end
            allstim = unique([unique(cl) unique(cl2)]);
            allstim2 = cellfun(@(x) fliplr(x),allstim,'uni',0);
            for istim = 1:length(Stims)
                idx = strcmp(Stims(istim),allstim2);
                if any(idx); Stims(istim)= allstim(idx);end
            end
            tr = all(arrayfun(@(x) any(strcmp(x,Stims)),cl,'uni',1));
            ts = all(arrayfun(@(x) any(strcmp(x,Stims)),cl2,'uni',1));
            if ~tr; disp 'Training stimuli not shown!'; return; end
            if ~ts; disp 'Testing stimuli not shown!'; return; end
                
            if ~isempty(test_groups)
                train_sz = cell2mat(cellfun(@size,train_groups,'uni',0));
                test_sz = cell2mat(cellfun(@size,test_groups,'uni',0));
                assert(all(train_sz==test_sz),'Test groups must match train groups')
            end
            
            % run the decoding
            P = cell(length(train_groups),length(train_groups{1})); P_shfl = P;
            for iGroup = 1:length(train_groups)
                train_data = [];test_data = [];
                for iClass = 1:length(train_groups{iGroup})
                    tgroup = train_groups{iGroup}{iClass};
                    stim_idx = any(strcmp(...
                        repmat(tgroup,size(Stims,1),1),...
                        repmat(Stims,1,size(tgroup,2)))',1);
                    train_data{iClass} = cell2mat(Traces(stim_idx));
                    if ~isempty(test_groups)
                        tgroup = test_groups{iGroup}{iClass};
                        stim_idx = any(strcmp(...
                            repmat(tgroup,size(Stims,1),1),...
                            repmat(Stims,1,size(tgroup,2)))',1);
                        test_data{iClass} = cell2mat(Traces(stim_idx));
                    end
                end
                [P(iGroup,:), P_shfl(iGroup,:)]= ...
                    decodeMulti(self,train_data,test_data,k_fold,shuffle,decoder,repetitions);
            end
            
            % insert
            key.p = P;
            key.p_shuffle = P_shfl;
            key.train_groups = train_groups;
            key.test_groups = test_groups;
            self.insert(key)
        end
    end
    
    methods
        
        function [train_groups, test_groups] = getGroups(obj,Stims,train_set,test_set)
            
            if isempty(train_set) % take all pairwise combinations of stimuli
                Stims = combnk(Stims,2);
                for igroup = 1:size(Stims,1)
                    for istim = 1:size(Stims,2)
                        train_groups{igroup}{istim} = Stims(igroup,istim);
                    end
                end
                test_groups = [];
            else
                train_groups = splitGroups(train_set);
                test_groups = splitGroups(test_set);
            end
            
            function groups = splitGroups(set)
                % output {group}{class}{obj}
                group_ids = regexp(set,'([^:\{\}$]*)','match');
                groups = [];
                for igroup = 1:length(group_ids)
                    classes = strsplit(group_ids{igroup},';');
                    for iclass = 1:length(classes)
                        groups{igroup}{iclass} = strsplit(classes{iclass},'.');
                    end
                end
            end
        end
        
        function [PP, RR] = decodeMulti(obj,Data,test_Data,k_fold,shuffle,decoder,repetitions)
            % performs a svm classification
            % data: {classes}[cells trials]
            % output: {classes}[reps trials]
            
            if nargin<5; shuffle = 0;end
            if nargin<4; k_fold = 10;end
            if nargin<3; test_Data = [];end
            
            % linearize
            Data = Data(:)';test_Data = test_Data(:)';
            
            PP = cell(repetitions,1); RR = PP;
            
            for irep = 1:repetitions
                % initialize
                groups = []; test_groups = []; train_idx = []; test_idx = [];
                s = RandStream('mt19937ar','Seed','shuffle');
                group_num = length(Data);
                
                % equalize by undersampling shorter class & randomize trial sequence
                msz = min(cellfun(@(x) size(x,2),Data));
                data = cellfun(@(x) x(:,randperm(s,size(x,2),msz)),Data,'uni',0);
                
                % use data as test_data if not provided
                if isempty(test_Data)
                    test_data = data;
                    s.reset;
                    data_idx = cellfun(@(x) randperm(s,size(x,2),msz),Data,'uni',0);% create bin index
                else % randomize trials
                    s = RandStream('mt19937ar','Seed','shuffle');
                    test_data = cellfun(@(x) x(:,randperm(s,size(x,2))),test_Data,'uni',0);
                    s.reset;
                    data_idx = cellfun(@(x) randperm(s,size(x,2)),test_Data,'uni',0);
                end
                
                % make group identities & build indexes
                if k_fold<2;bins = msz;else;bins = k_fold;end
                bin_sz = max([floor(msz/bins) 1]);
                
                for iclass = 1:group_num
                    % make group identities
                    groups{iclass} = ones(1,size(data{iclass},2)) * iclass;
                    test_groups{iclass} = ones(1,size(test_data{iclass},2)) * iclass;
                    
                    % buld index
                    test_bin_sz = max([floor(size(test_data{iclass},2)/bins) 1]);
                    for ibin = 1:bins
                        train_idx{iclass}(1 + (ibin-1)*bin_sz:bin_sz*ibin) = ibin;
                        test_idx{iclass}(1 + (ibin-1)*test_bin_sz:test_bin_sz*ibin) = ibin;
                    end
                end
                
                % combine classes in one vector
                data = cell2mat(data);
                groups = cell2mat(groups);
                test_data = cell2mat(test_data);
                test_groups = cell2mat(test_groups);
                train_idx = cell2mat(train_idx);
                test_idx = cell2mat(test_idx);
                data_idx = cell2mat(data_idx);
                
                % create shuffled testing trials
                test_sz = size(test_data,2);
                test_shfl_idx = 1:size(test_groups,2);
                for ishuffle = 1:shuffle
                    test_shfl_idx = test_shfl_idx(randperm(test_sz));
                end
                test_shfl_groups = test_groups(test_shfl_idx);
                data_shfl_idx = data_idx(test_shfl_idx);
                
                % classify
                P = cellfun(@(x) nan(1,size(x,2)),Data,'uni',0);R = P;
                for ibin = 1:bins
                    idx = train_idx ~= ibin;
                    tidx = test_idx == ibin;
                    DEC = feval(decoder,data(:,idx)', groups(idx)');
                    pre = predict(DEC,test_data(:,tidx)');
                    p =  (pre == test_groups(tidx)');
                    r =  (pre == test_shfl_groups(tidx)');
                    for igroup = 1:group_num
                        P{igroup}(data_idx(tidx & test_groups==igroup)) = p(test_groups(tidx)==igroup);
                        R{igroup}(data_shfl_idx(tidx & test_shfl_groups==igroup)) = r(test_shfl_groups(tidx)==igroup);
                    end
                end
                PP{irep} = P;
                RR{irep} = R;
            end
            
            % convert {reps}{obj}[1 trials] to {obj}[reps trials]
            PP = cellfun(@cell2mat,mat2cell(permute(reshape([PP{:}],...
                length(Data),repetitions),[2 1]),repetitions,ones(1,length(Data))),'uni',0);
            RR = cellfun(@cell2mat,mat2cell(permute(reshape([RR{:}],...
                length(Data),repetitions),[2 1]),repetitions,ones(1,length(Data))),'uni',0);
        end
        
    end
end
