%{
# Object decoding
-> fuse.ScanDone
-> anatomy.Area
-> obj.DecodeOpt
---
-> stimulus.Sync
p                     : longblob                      # performance
p_shuffle             : longblob                      # chance performance
train_groups          : longblob                      # train group info
test_groups           : longblob                      # test group info
trial_info            : longblob                      # trial info [index, clip #, movie, bin #]
%}

classdef Decode < dj.Computed
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.DecodeOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"'))))
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            [decoder,k_fold,shuffle,train_set,test_set,repetitions,select_method] = ...
                fetch1(obj.DecodeOpt & key,...
                'decoder','k_fold','shuffle','train_set','test_set','repetitions','select_method');
            
            [Traces, Stims, StimInfo] = getData(self,key); % [Cells, Obj, Trials]
            [train_groups,test_groups] = getGroups(self,Stims,train_set,test_set);
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
                    info.bins{iGroup,iClass} = cell2mat(StimInfo.bins(stim_idx));
                    info.trials{iGroup,iClass} = cell2mat(StimInfo.trials(stim_idx));
                    info.clips{iGroup,iClass} = cell2mat(StimInfo.clips(stim_idx));
                    info.names{iGroup,iClass} = [StimInfo.names{stim_idx}];
                end
                
                [P(iGroup,:), P_shfl(iGroup,:)]= ...
                    decodeMulti(self,train_data,test_data,k_fold,shuffle,decoder,repetitions,select_method);
            end
            
            % insert
            key.p = P;
            key.p_shuffle = P_shfl;
            key.train_groups = train_groups;
            key.test_groups = test_groups;
            key.trial_info = info;
            self.insert(key)
        end
    end
    
    methods
        function [Data, Stims, info] = getData(self,key,bin)
            
            if nargin<3
                bin = fetch1(obj.DecodeOpt & key, 'binsize');
            end
            
            % get traces
            [Traces, caTimes] = getAdjustedSpikes(fuse.ActivityTrace & key,'soma'); % [time cells]
            
            % get rid of nans
            notnanidx = ~isnan(mean(Traces,2)); % faster than all
            Traces = Traces(notnanidx,:);
            caTimes = caTimes(notnanidx);
            
            % interpolate over time
            X = @(t) interp1(caTimes-caTimes(1), Traces, t, 'linear', 'extrap');  % traces indexed by time
            
            % fetch stuff
            [flip_times, trial_idxs] = fetchn(stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"')) & key,'flip_times','trial_idx');
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            trial_idxs = trial_idxs(tidx);
            flip_times = cell2mat(flip_times(tidx));
            Stims = unique(fetchn(stimulus.Clip &  (stimulus.Trial & key),'movie_name'));
            
            % subsample traces
            fps = 1/median(diff(flip_times(1,:)));
            d = max(1,round(bin/1000*fps));
            Traces = convn(permute(X(flip_times - caTimes(1)),[2 3 1]),ones(d,1)/d,'same');
            Traces = permute(Traces(1:d:end,:,:),[2 1 3]); % in [cells bins trials]
            
            % split for unique stimuli
            for istim = 1:length(Stims)
                [s_trials,s_clips,s_names] = fetchn(stimulus.Trial * stimulus.Clip & ...
                    sprintf('movie_name = "%s"',Stims{istim}) & key, 'trial_idx','clip_number','movie_name');
                [tr_idx, b]= ismember(trial_idxs,s_trials);
                st_idx = b(b>0);
                dat = Traces(:,:,tr_idx);
                info.bins{istim} = reshape(repmat(1:size(dat,2),size(dat,3),1)',[],1);
                info.trials{istim} = reshape(repmat(s_trials(st_idx),1,size(dat,2))',[],1);
                info.clips{istim} = reshape(repmat(s_clips(st_idx),1,size(dat,2))',[],1);
                info.names{istim} = reshape(repmat(s_names(st_idx),1,size(dat,2))',[],1);
                Data{istim} = reshape(dat,size(Traces,1),[]);
            end
        end
        
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
                        groups{igroup}{iclass} = strsplit(classes{iclass},',');
                    end
                end
            end
        end
        
        function [PP, RR] = decodeMulti(obj,Data,test_Data,k_fold,shuffle,decoder,repetitions,select_method)
            % performs a svm classification
            % data: {classes}[cells trials]
            % output: {classes}[reps trials]
            
            if nargin<5; shuffle = 0;end
            if nargin<4; k_fold = 10;end
            if nargin<3; test_Data = [];end
            
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
                bin_sz = floor(msz/bins);
                
                for iclass = 1:group_num
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
                
                % combine classes in one vector
                data = cell2mat(data);
                groups = cell2mat(groups);
                test_data = cell2mat(test_data);
                test_groups = cell2mat(test_groups);
                train_idx = cell2mat(train_idx);
                test_idx = cell2mat(test_idx);
                data_idx = cell2mat(data_idx);
                
                % make nan zeros
                data(isnan(data)) = prctile(data(:),1);
                
                % create shuffled testing trials
                test_sz = size(test_data,2);
                test_shfl_idx = 1:size(test_groups,2);
                for ishuffle = 1:shuffle
                    test_shfl_idx = test_shfl_idx(randperm(test_sz));
                end
                test_shfl_groups = test_groups(test_shfl_idx);
                data_shfl_idx = data_idx(test_shfl_idx);
                
                % get cell index
                cell_idx = randperm(size(data,1));
                if strcmp(select_method,'subsample')
                    cell_num = 1:2:size(data,1);
                else
                    cell_num = size(data,1);
                end
                
                % classify
                P = cellfun(@(x) nan(length(cell_num),size(x,2)),Data,'uni',0);R = P;
                for icell = 1:length(cell_num)
                    icelln = cell_num(icell);
                    for ibin = 1:bins
                        idx = train_idx ~= ibin;
                        tidx = test_idx == ibin;
                        DEC = feval(decoder,data(cell_idx(1:icelln),idx)', groups(idx)','learner','svm',...
                            'regularization','lasso','solver','sparsa');
                        pre = predict(DEC,test_data(cell_idx(1:icelln),tidx)');
                        p =  (pre == test_groups(tidx)');
                        r =  (pre == test_shfl_groups(tidx)');
                        for igroup = 1:group_num
                            P{igroup}(icell,data_idx(tidx & test_groups==igroup)) = p(test_groups(tidx)==igroup);
                            R{igroup}(icell,data_shfl_idx(tidx & test_shfl_groups==igroup)) = r(test_shfl_groups(tidx)==igroup);
                        end
                    end
                end
                PP{irep} = P;
                RR{irep} = R;
            end
            
            % convert {reps}{obj}[cells trials] to {obj}[reps trials cells]
            PP = cellfun(@cell2mat,mat2cell(cellfun(@(x) permute(x,[3 2 1]),permute(reshape([PP{:}],...
                length(Data),repetitions),[2 1]),'uni',0),repetitions,ones(1,length(Data))),'uni',0);
            RR = cellfun(@cell2mat,mat2cell(cellfun(@(x) permute(x,[3 2 1]),permute(reshape([RR{:}],...
                length(Data),repetitions),[2 1]),'uni',0),repetitions,ones(1,length(Data))),'uni',0);
            
        end
        
        function plotMasks(self,norm)
            
            % get data
            method = fetch1(obj.DecodeOpt & self,'decoder');
            [mi, area] = fetchn(self,'p','brain_area');
            areas = unique(area);
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                MI{iarea} = cellfun(@(x) nanmean(reshape(cellfun(@(xx) nanmean(xx(:)),x),[],1)), mi(strcmp(area,areas(iarea))));
            end
            
            % plot
            f = figure;
            colors = parula(30);
            plotMask(anatomy.Masks)
            colormap parula
            c = colorbar;
            name = 'Classification performance (%)';
            ylabel(c,name,'Rotation',-90,'VerticalAlignment','baseline')
            
            mxMI = max(cellfun(@nanmean,MI));
            if nargin>1
                mx = mxMI;
            else
                mx =1;
            end
            for iarea = 1:length(areas)
                mi = nanmean(MI{iarea});
                idx = double(uint8(floor(((mi-0.5)/(mx - 0.5))*0.99*size(colors,1)))+1);
                plotMask(anatomy.Masks & ['brain_area="' areas{iarea} '"'],colors(idx,:),length(MI{iarea}))
            end
            set(c,'ytick',linspace(0,1,5),'yticklabel',roundall(linspace(0.5,mx,5),0.01))
        end
    end
end
