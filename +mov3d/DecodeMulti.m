%{
mov3d.DecodeMulti (computed) # calcium trace
-> preprocess.Sync
-> mov3d.DecodeMultiOpt
-> preprocess.SpikeMethod
-> preprocess.Method
---
p                     : longblob                      # performance
p_shuffle             : longblob                      # chance performance
train_groups          : longblob                      # train group info
test_groups           : longblob                      # test group info
%}

classdef DecodeMulti < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  = (experiment.Scan  ...
            * (preprocess.Spikes))...
            * (mov3d.DecodeMultiOpt & 'process = "yes"') ...
            * (preprocess.Sync  & (vis.MovieClipCond * vis.Trial & ...
            (vis.Movie & 'movie_class="multiobjects"') & ...
            'trial_idx between first_trial and last_trial'))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            [decoder,k_fold,shuffle,train_set,test_set,repetitions] = ...
                fetch1(mov3d.DecodeMultiOpt & key,...
                'decoder','k_fold','shuffle','train_set','test_set','repetitions');
            
            [Traces, Stims] = getData(self,key); % [Cells, Obj, Trials]
            %stim_obj = cellfun(@cell2mat,regexp(Stims,'(?<=(o))\d+','match'),'uni',0);
            train_groups = getGroups(self,train_set);
            test_groups = getGroups(self,test_set);
            if ~isempty(test_groups)
                train_sz = cell2mat(cellfun(@size,train_groups,'uni',0));
                test_sz = cell2mat(cellfun(@size,test_groups,'uni',0));
                assert(all(train_sz==test_sz),'Test groups must match train groups')
            end
            
            % run the decoding
            P = nan(length(train_groups),k_fold,repetitions); P_shfl = P;
            for iGroup = 1:length(train_groups)
                train_data = [];test_data = [];
                for iClass = 1:length(train_groups{iGroup})
                    tgroup = train_groups{iGroup}{iClass};
                    stim_idx = any(strcmp(...
                        repmat(tgroup,size(Stims,1),1),...
                        repmat(Stims,1,size(tgroup,2)))');
                    train_data{iClass} = cell2mat(Traces(stim_idx));
                    if ~isempty(test_groups)
                        tgroup = test_groups{iGroup}{iClass};
                        stim_idx = any(strcmp(...
                            repmat(tgroup,size(Stims,1),1),...
                            repmat(Stims,1,size(tgroup,2)))');
                        test_data{iClass} = cell2mat(Traces(stim_idx));
                    end
                end
                [P(iGroup,:,:), P_shfl(iGroup,:,:)]= ...
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
        function [Data, Stims] = getData(obj,key,ibin)
            
            bin = fetch1(mov3d.DecodeMultiOpt & key, 'binsize');
            if nargin>2;bin = ibin;end
            
            % get traces
            [Traces, caTimes] = pipetools.getAdjustedSpikes(key);
            xm = min([length(caTimes) length(Traces)]);
            X = @(t) interp1(caTimes(1:xm)-caTimes(1), Traces(1:xm,:), t, 'linear', nan);  % traces indexed by time
            
            % fetch stuff
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            Stims = unique(fetchn(vis.MovieClipCond & trials,'movie_name'));
            [flip_times, trial_idxs] = (fetchn(trials * vis.MovieClipCond,'flip_times','trial_idx'));
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            trial_idxs = trial_idxs(tidx);
            flip_times = cell2mat(flip_times(tidx));
            
            % subsample traces
            fps = 1/median(diff(flip_times(1,:)));
            d = max(1,round(bin/1000*fps));
            traces = convn(permute(X(flip_times - caTimes(1)),[2 3 1]),ones(d,1)/d,'same');
            traces = traces(1:d:end,:,:);
            
            % split for unique stimuli
            for istim = 1:length(Stims)
                k.movie_name = Stims{istim};
                stim_trials = fetchn(trials*vis.MovieClipCond & k,'trial_idx');
                Data{istim} = reshape(permute(traces(:,:,ismember(trial_idxs,stim_trials)),[2 1 3]),size(traces,2),[]);
            end
        end
        
        function groups = getGroups(obj,set)
            % output {group}{class}{obj}
            group_ids = regexp(set,'([^:\{\}$]*)','match');
            groups = [];
            for igroup = 1:length(group_ids)
                classes = strsplit(group_ids{igroup},';');
                for iclass = 1:length(classes)
                    groups{igroup}{iclass} = strsplit(classes{iclass},',');
                end
            end
        end
        
        function [PP, RR] = decodeMulti(obj,Data,test_Data,k_fold,shuffle,decoder,repetitions)
            % performs a svm classification
            % data: {classes}[cells trials]

            if nargin<5; shuffle = 0;end
            if nargin<4; k_fold = 10;end
            if nargin<3; test_Data = [];end
            
            parfor irep = 1:repetitions
                % initialize
                groups = []; test_groups = []; train_idx = []; test_idx = [];

                % equalize by undersampling shorter class & randomize trial sequence
                msz = min(cellfun(@(x) size(x,2),Data));
                data = cellfun(@(x) x(:,randperm(size(x,2),msz)),Data,'uni',0);

                % use data as test_data if not provided
                if isempty(test_Data)
                    test_data = data;
                else % randomize trials
                    test_data = cellfun(@(x) x(:,randperm(size(x,2))),test_Data,'uni',0);
                end

                % make group identities & build indexes
                if k_fold<2;bins = msz;else;bins = k_fold;end
                bin_sz = floor(msz/bins);

                for iclass = 1:length(data)
                    % make group identities
                    groups{iclass} = ones(1,size(data{iclass},2)) * iclass;
                    test_groups{iclass} = ones(1,size(test_data{iclass},2)) * iclass;

                    % buld index
                    test_bin_sz = floor(size(test_data{iclass},2)/bins);
                    for ibin = 1:bins
                       train_idx{iclass}(1 + (ibin-1)*bin_sz:bin_sz*ibin) = ibin;
                       test_idx{iclass}(1 + (ibin-1)*test_bin_sz:test_bin_sz*ibin) = ibin;
                    end
                end

                % make data vectors
                data = cell2mat(data);
                groups = cell2mat(groups);
                test_data = cell2mat(test_data);
                test_groups = cell2mat(test_groups);
                train_idx = cell2mat(train_idx);
                test_idx = cell2mat(test_idx);
                test_sz = size(test_data,2);

                % create shuffled testing trials
                test_shfl_groups = test_groups;
                for ishuffle = 1:shuffle
                    test_shfl_groups = test_shfl_groups(randperm(test_sz));
                end

                % classify
                P = nan(bins,1);R = P;
                for ibin = 1:bins
                    idx = train_idx ~= ibin;
                    tidx = test_idx == ibin;
                    DEC = feval(decoder,data(:,idx)', groups(idx)'); 
                    pre = predict(DEC,test_data(:,tidx)');
                    P(ibin) =  mean(pre == test_groups(tidx)');
                    R(ibin) =  mean(pre == test_shfl_groups(tidx)');
                end
                PP(:,irep) = P;
                RR(:,irep) = R;
            end
        end
        
        function plotMasks(obj,varargin)
            params.normalize = 0;
            params.colors = 30;
            params.steps =11;
            params = getParams(params,varargin);
            
            figure
            colors = parula(params.colors);
            plotMask(map.Masks)
            colormap parula
            c = colorbar;
            ylabel(c,'Performance','Rotation',-90,'VerticalAlignment','baseline')
            areas =  fetchn(map.Area,'area');
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                keys = fetch(obj & (experiment.Scan & ['brain_area="' areas{iarea} '"']));
                if isempty(keys);continue;end 
                for ikey = 1:length(keys)
                    tuple = keys(ikey);
                    MI{iarea}(ikey) = mean(reshape(fetch1(obj & tuple,'p'),[],1));
                end
            end
            
            if params.normalize
                mx = max(cell2mat(cellfun(@max,MI,'uni',0)));
                mn = min(cell2mat(cellfun(@min,MI,'uni',0)));
                params.steps = 5;
            else
                mn = 0.5;
                mx = 1;
            end
            labels = round(linspace(mn,mx,params.steps)*100)/100;
            color_idx = @(x) round(min([max([(x-mn)/(mx-mn) 0]) 1])*(params.colors-1))+1;

            for iarea = 1:length(areas)
                mi = mean(MI{iarea});
                if isnan(mi);continue;end
                idx = color_idx(mi);
                plotMask(map.Masks & ['area="' areas{iarea} '"'],colors(idx,:),length(MI{iarea}))
            end

            set(c,'ytick',linspace(0,1,length(labels)),'yticklabel',labels)

            
        end
    end
end
