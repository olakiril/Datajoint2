%{
# Object decoding
-> fuse.ScanDone
-> obj.DecodeOpt
brain_area_1 : varchar(12)
brain_area_2 : varchar(12)
---
-> stimulus.Sync
p                     : longblob                      # performance
p_shuffle             : longblob                      # chance performance
train_groups          : longblob                      # train group info
test_groups           : longblob                      # test group info
trial_info            : longblob                      # trial info [index, clip #, movie, bin #]
score                 : longblob                      # svm score (distance to boundary)
classifier            : longblob                      # svm score (distance to boundary)
%}

classdef DecMulti < dj.Computed
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        keySource  = fuse.ScanDone...
            * (obj.DecodeOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  ...
            (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"'))))
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % get train & test groups
            [train_set,test_set, neurons] = fetch1(obj.DecodeOpt & key,'train_set','test_set','neurons');
            
            areas = fetchn(fuse.ScanDone * anatomy.Area & anatomy.AreaMembership & key,'brain_area');
            
            All_traces = cell(length(areas),1);
            All_unit_ids = cell(length(areas),1);
            for iarea = 1:length(areas)
                % get Datad
                area_key = key;
                area_key.brain_area = areas{iarea};
                [All_traces{iarea}, Stims, StimInfo, All_unit_ids{iarea}, BehTraces] = getData(obj.Dec,area_key); % [Cells, Obj, Trials]
            end
            
            [train_groups,test_groups] = getGroups(obj.Dec,Stims,train_set,test_set);
            if ~isempty(test_groups)
                train_sz = cell2mat(cellfun(@size,train_groups,'uni',0));
                test_sz = cell2mat(cellfun(@size,test_groups,'uni',0));
                assert(all(train_sz==test_sz),'Test groups must match train groups')
            end
            
            % merge selection of cells
            sz = cellfun(@(x) size(x{1},1),All_traces);
            idx = sz>neurons/2; % select only areas that have half the required neurons
            assert(sum(idx)>1,'Not enough neurons for comparisson!')
            All_traces = All_traces(idx);
            All_unit_ids = All_unit_ids(idx);
            areas = areas(idx);
            area_combs = combnk(1:sum(idx),2);
            
            for icomb = 1:size(area_combs,1)
                mn_sz = min(sz(area_combs(icomb,:)));
                cell_idx1 = randperm(sz(area_combs(icomb,1)),mn_sz);
                cell_idx2 = randperm(sz(area_combs(icomb,2)),mn_sz);
                Traces = cellfun(@(x,y) [x(cell_idx1,:);y(cell_idx2,:)],All_traces{area_combs(icomb,1)},All_traces{area_combs(icomb,2)},'UniformOutput',false);
                Unit_ids = [All_unit_ids{area_combs(icomb,1)}(cell_idx1),All_unit_ids{area_combs(icomb,2)}(cell_idx2)];
                
                %Loop through each group of comparissons
                P = cell(length(train_groups),length(train_groups{1})); P_shfl = P;score = P;unit_idx = {};classifier = {};
                for iGroup = 1:length(train_groups)
                    fprintf('\n Group# %d/%d ',iGroup,length(train_groups));
                    
                    % assign training & testing data and the stimulus information for each
                    train_data = [];test_data = [];class_idx = true(length(train_groups{iGroup}),1);
                    beh_data = [];
                    for iClass = 1:length(train_groups{iGroup})
                        
                        % Training Data
                        tgroup = train_groups{iGroup}{iClass};
                        stim_idx = any(strcmp(...
                            repmat(tgroup,size(Stims,1),1),...
                            repmat(Stims,1,size(tgroup,2)))',1);
                        if all(stim_idx==0);class_idx(iClass) = false;continue;end
                        train_data{iClass} = cell2mat(Traces(stim_idx));
                        beh_data{iClass} = cell2mat(BehTraces(stim_idx));
                        train_info.bins{iGroup,iClass} = cell2mat(reshape(StimInfo.bins(stim_idx),[],1));
                        train_info.trials{iGroup,iClass} = cell2mat(reshape(StimInfo.trials(stim_idx),[],1));
                        train_info.clips{iGroup,iClass} = cell2mat(reshape(StimInfo.clips(stim_idx),[],1));
                        names = cellfun(@(x) reshape(x,1,[]),StimInfo.names,'uni',0);
                        train_info.names{iGroup,iClass} = [names{stim_idx}];
                        
                        % Testing Data
                        if ~isempty(test_groups)
                            tgroup = test_groups{iGroup}{iClass};
                            stim_idx = any(strcmp(...
                                repmat(tgroup,size(Stims,1),1),...
                                repmat(Stims,1,size(tgroup,2)))',1);
                            test_data{iClass} = cell2mat(Traces(stim_idx));
                            test_info.bins{iGroup,iClass} = cell2mat(reshape(StimInfo.bins(stim_idx),[],1));
                            test_info.trials{iGroup,iClass} = cell2mat(reshape(StimInfo.trials(stim_idx),[],1));
                            test_info.clips{iGroup,iClass} = cell2mat(reshape(StimInfo.clips(stim_idx),[],1));
                            names = cellfun(@(x) reshape(x,1,[]),StimInfo.names,'uni',0);
                            test_info.names{iGroup,iClass} = [names{stim_idx}];
                        else, test_info = train_info;
                        end
                    end
                    if isempty(train_data); error('No data selected!');end
                    
                    % run the decoding
                    [P(iGroup,class_idx), P_shfl(iGroup,class_idx), unit_idx(iGroup,:), score(iGroup,class_idx), classifier(iGroup,:)]= ...
                        decodeMulti(obj.Dec, train_data, test_data, key, Unit_ids,...
                        structfun(@(x) x(iGroup,:),train_info,'uni',0),structfun(@(x) x(iGroup,:),test_info,'uni',0), beh_data);
                end
                
                % find unit ids from randomization indexes
                test_info.units = cellfun(@(x) cellfun(@(xx) Unit_ids(xx),x,'uni',0),unit_idx,'uni',0);
                
                % insert
                key.p = P;
                key.p_shuffle = P_shfl;
                key.train_groups = train_groups;
                key.test_groups = test_groups;
                key.trial_info = test_info;
                key.score = score;
                key.classifier = classifier;
                key.p_shuffle = {};
                key.brain_area_1 = areas{area_combs(icomb,1)};
                key.brain_area_2 = areas{area_combs(icomb,2)};
                for i = 1:length(key.score)
                    key.score{i} = int16(key.score{i}*10000);
                    key.p{i} = int8(key.p{i});
                end
                self.insert(key)
            end
        end
    end
end