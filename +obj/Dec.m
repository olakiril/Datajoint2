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
score                 : longblob                      # svm score (distance to boundary)
classifier            : longblob                      # svm score (distance to boundary)
%}

classdef Dec < dj.Computed
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.DecodeOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  ...
            (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"'))))
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % get train & test groups
            [train_set,test_set] = fetch1(obj.DecodeOpt & key,'train_set','test_set');
            
            % get Data
            [Traces, Stims, StimInfo, Unit_ids] = getData(self,key); % [Cells, Obj, Trials]
            [train_groups,test_groups] = getGroups(self,Stims,train_set,test_set);
            if ~isempty(test_groups)
                train_sz = cell2mat(cellfun(@size,train_groups,'uni',0));
                test_sz = cell2mat(cellfun(@size,test_groups,'uni',0));
                assert(all(train_sz==test_sz),'Test groups must match train groups')
            end
            
            %Loop through each group of comparissons
            P = cell(length(train_groups),length(train_groups{1})); P_shfl = P;score = P;unit_idx = {};classifier = {};
            for iGroup = 1:length(train_groups)
                fprintf('Group# %d/%d ',iGroup,length(train_groups));
                
                % assign training & testing data and the stimulus information for each
                train_data = [];test_data = [];class_idx = true(length(train_groups{iGroup}),1);
                for iClass = 1:length(train_groups{iGroup})
                    
                    % Training Data
                    tgroup = train_groups{iGroup}{iClass};
                    stim_idx = any(strcmp(...
                        repmat(tgroup,size(Stims,1),1),...
                        repmat(Stims,1,size(tgroup,2)))',1);
                    if all(stim_idx==0);class_idx(iClass) = false;continue;end
                    train_data{iClass} = cell2mat(Traces(stim_idx));
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
                    decodeMulti(self, train_data, test_data, key, Unit_ids,...
                    structfun(@(x) x(iGroup,:),train_info,'uni',0),structfun(@(x) x(iGroup,:),test_info,'uni',0));
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
            self.insert(key)
        end
    end
    
    methods
        function [PP, RR, Cells, SC, DC] = decodeMulti(self,Data,test_Data, key, unit_ids, train_info, test_info)
            % performs a svm classification
            % data: {classes}[cells trials]
            % output: {classes}[reps trials]
            
            % get decoder parameters
            [decoder,k_fold,shuffle,repetitions,select_method, ...
                dec_params, neurons,fold_selection, binsize, dat_norm] = ...
                fetch1(obj.DecodeOpt & key,...
                'decoder','k_fold','shuffle','repetitions','select_method',...
                'dec_params','neurons','fold_selection','binsize','normalize');
            
            
            % define decoder function
            if ~isempty(dec_params);dec_params = [',' dec_params];end
            decoder_func = eval(sprintf('@(X,IDs) %s(X, IDs%s)',decoder,dec_params));
            
            % get test data
            if isempty(test_Data)
                test_Data = Data;
            end
            
            PP = cell(repetitions,1); RR = PP;Cells = [];SC = PP;txt = '';
            for irep = 1:repetitions
                rseed = RandStream('mt19937ar','Seed',irep);
                
                % initialize
                groups = []; test_groups = []; train_idx = []; test_idx = [];
                group_num = length(Data);
                
                % equalize by undersampling shorter class & randomize trial sequence
                msz = min(cellfun(@(x) size(x,2),Data)); % calculate minimum class length
                
                % calculate fold bin size and recompute minimum length of data
                if k_fold==1;bins = msz;elseif k_fold==0;bins=1;else, bins = k_fold;end
                bin_sz = floor(msz/bins);
                msz = bin_sz*bins;
                
                % undersample data
                train_data_idx = cellfun(@(x) randperm(rseed,size(x,2),msz),Data,'uni',0);% create bin index
                
                % equalize by undersampling shorter class & randomize trial sequence
                msz = min(cellfun(@(x) size(x,2),test_Data)); % calculate minimum class length
                test_bin_sz = floor(msz/bins); % calculate fold bin size and recompute minimum length of data
                msz = test_bin_sz*bins;
                rseed.reset;
                test_data_idx = cellfun(@(x) randperm(rseed,size(x,2),msz),test_Data,'uni',0);
                
                % make group identities & build indexes
                for iclass = 1:group_num
                    % make group identities
                    groups{iclass} = ones(1,size(train_data_idx{iclass},2)) * iclass;
                    test_groups{iclass} = ones(1,size(test_data_idx{iclass},2)) * iclass;
                    
                    switch fold_selection
                        case 'random'
                            % buld index
                            for ibin = 1:bins
                                train_idx{iclass}(1 + (ibin-1)*bin_sz      :      bin_sz*ibin) = ibin;
                                test_idx{iclass} (1 + (ibin-1)*test_bin_sz : test_bin_sz*ibin) = ibin;
                            end
                        case {'x','y','frame_id','scale','rotation','tilt','light1_xloc','light1_yloc','light3_yloc','light2_ene'}
                            
                            keys = cell2struct([num2cell(train_info.clips{iclass}(train_data_idx{iclass}),2),...
                                train_info.names{iclass}(train_data_idx{iclass})']',{'clip_number','movie_name'},1);
                            param = getParam(stimulus.MovieParams, keys, fold_selection, ...
                                train_info.bins{iclass}(train_data_idx{iclass})*abs(binsize)/1000 - abs(binsize)/2/1000);
                            [~, sort_idx] = histc(param,[-inf; quantile(param,bins-1)'; inf]);
                            train_idx{iclass} = sort_idx(:)';
                            keys = cell2struct([num2cell(test_info.clips{iclass}(test_data_idx{iclass}),2),...
                                test_info.names{iclass}(test_data_idx{iclass})']',{'clip_number','movie_name'},1);
                            param = getParam(stimulus.MovieParams, keys, fold_selection, ...
                                test_info.bins{iclass}(test_data_idx{iclass})*abs(binsize)/1000 - abs(binsize)/2/1000);
                            [~, sort_idx] = histc(param,[-inf; quantile(param,bins-1)'; inf]);
                            test_idx{iclass} = sort_idx(:)';
                        case 'time'
                              % buld index
                            for ibin = 1:bins
                                train_idx{iclass}(1 + (ibin-1)*bin_sz      :      bin_sz*ibin) = ibin;
                                test_idx{iclass} (1 + (ibin-1)*test_bin_sz : test_bin_sz*ibin) = ibin;
                            end
                            
                            bns = test_info.bins{iclass};
                            tls = test_info.trials{iclass};
                            un_trl = unique(tls);
                            for trial = un_trl(:)'
                                idx = tls==trial;
                                [bi,sbi] = sort(bns(idx),'ASC');
                                dat = test_Data{iclass}(:,idx);
                                dat = dat(:,sbi);
                                dat(:,2:max(bi)) = dat(:,1:max(bi)-1);
                                dat(:,1) = nan;
                                test_Data{iclass}(:,idx) = dat(:,sbi);
                            end
                          
                        otherwise
                            error('Fold selection method not implemented!')
                    end
                end
                
                % combine classes in one vector
                data = cell2mat(cellfun(@(x,y) x(:,y),Data,train_data_idx,'uni',0));
                groups = cell2mat(groups);
                test_data = cell2mat(cellfun(@(x,y) x(:,y),test_Data,test_data_idx,'uni',0));
                test_groups = cell2mat(test_groups);
                train_idx = cell2mat(train_idx);
                test_idx = cell2mat(test_idx);
                test_data_idx = cell2mat(test_data_idx);
                
                % normalize data
                if dat_norm
                    mdata = mean(data,2);
                    sdata = std(data,[],2);
                else
                    mdata = zeros(size(data,1),1);
                    sdata = ones(size(data,1),1);
                end
                
                data = bsxfun(@rdivide, bsxfun(@minus,data,mdata),sdata);
                test_data = bsxfun(@rdivide, bsxfun(@minus,test_data,mdata),sdata);
                
                % make nan zeros
                data(isnan(data)) = prctile(data(:),1);
                
                % create shuffled testing trials
                test_sz = size(test_data,2);
                test_shfl_idx = 1:size(test_groups,2);
                for ishuffle = 1:shuffle
                    test_shfl_idx = test_shfl_idx(randperm(rseed,test_sz));
                end
                test_shfl_groups = test_groups(test_shfl_idx);
                data_shfl_idx = test_data_idx(test_shfl_idx);
                
                % get cell index
                cell_idx = randperm(rseed,size(Data{1},1));
                switch select_method
                    case 'subsample'
                        cell_num = tril(true(numel(cell_idx)),0);
                        cell_num = cell_num([1 2.^(1:log2(numel(cell_idx)-1)) numel(cell_idx)],:);
                    case 'subsample2'
                        ncells = numel(cell_idx);
                        cell_num = tril(true(ncells),0);
                        cell_num = cell_num([floor(1.5.^(1:(log(ncells) / log(1.5)))) ncells],:);
                    case 'all'
                        cell_num = true(size(cell_idx));
                    case 'single'
                        cell_num = diag(true(size(cell_idx)));
                    case 'fixed'
                        if neurons>size(Data{1},1); fprintf('Recording only has %d neurons!',size(Data{1},1));return;end
                        cell_num = false(1,numel(cell_idx));
                        cell_num(1:neurons) = true;
                    case 'rf'
                        % get all rfs to compute center
                        [x,y] = getDegrees(tune.DotRFMap & (fuse.ScanDone & key) & 'p_value<0.05',1);
                        mx = mean(x);
                        my = mean(y);
                        [x,y, keys] = getDegrees(tune.DotRFMap & (fuse.ScanDone & key) & ...
                            'p_value<0.05' & (anatomy.AreaMembership & key),1);
                        idx = (mx-10)<x & x<(mx+10) & (my-10)<y & y<(my+10);
                        sel_units = [keys(idx).unit_id];
                        [~,unit_idx] = intersect(unit_ids(cell_idx),sel_units);
                        %indexes = [1 10:10:99 100:100:length(unit_idx) length(unit_idx)];
                        if neurons>size(Data{1},1); error('Recording only has %d neurons!',size(Data{1},1));end
                        indexes = neurons;
                        cell_num = false(length(indexes),numel(cell_idx));
                        for i = 1:length(indexes)
                            cell_num(i,unit_idx(1:indexes(i))) = true;
                        end
                    case 'reliable'
                        [r,keys] = fetchn(obj.RepeatsUnit & (anatomy.AreaMembership & key) & 'rep_opt = 2','r');
                        ids = [keys.unit_id];
                        un_ids = unique(ids);
                        rel = nan(length(un_ids),1);
                        for i = 1:length(un_ids)
                            rel(i) = nanmean((r(ids==un_ids(i))));
                        end
                        
                        % select nneurons most reliable
                        if neurons>size(Data{1},1); error('Recording only has %d neurons!',size(Data{1},1));end
                        [~,sort_idx] = sort(rel,'descend');
                        sel_units = un_ids(sort_idx);
                        indexes = neurons;
                        [~,unit_idx] = intersect(unit_ids(cell_idx),sel_units(1:indexes));
                        cell_num = false(length(indexes),numel(cell_idx));
                        for i = 1:length(indexes)
                            cell_num(i,unit_idx(1:indexes(i))) = true;
                        end
                    case 'reliablev2'
                        [r,keys] = fetchn(obj.RepeatsUnit & (anatomy.AreaMembership & key) & 'rep_opt = 2','r');
                        ids = [keys.unit_id];
                        un_ids = unique(ids);
                        rel = nan(length(un_ids),1);
                        for i = 1:length(un_ids)
                            rel(i) = nanmean((r(ids==un_ids(i))));
                        end
                        
                        % select nneurons with reliability above 0.2
                        if neurons>size(Data{1},1); error('Recording only has %d neurons!',size(Data{1},1));end
                        sel_units = un_ids(rel>0.2);
                        indexes = neurons;
                        [~,unit_idx] = intersect(unit_ids(cell_idx),sel_units);
                        cell_num = false(length(indexes),numel(cell_idx));
                        for i = 1:length(indexes)
                            cell_num(i,unit_idx(randperm(rseed,length(unit_idx),indexes))) = true;
                        end
                    otherwise
                        error('Cell selection method not supported')
                end
                
                % classify
                P = cellfun(@(x) nan(size(cell_num,1),size(x,2),'single'),test_Data,'uni',0);R = P;S = P;D = [];
                for icell = 1:size(cell_num,1) % For each cell permutation
                    for ibin = 1:bins          % For each fold
                        % print report
                        fprintf(repmat('\b',1,length(txt)));
                        txt= sprintf('Rep# %d/%d  Cell# %d/%d  Bin# %d/%d',irep,repetitions,icell,size(cell_num,1),ibin,bins);
                        fprintf('%s',txt);
                        
                        % select training/testing bin
                        if bins==1 % no crossvalidation case, testing on training dataset
                            idx = logical(train_idx);
                            tidx = logical(test_idx);
                        else
                            idx = train_idx ~= ibin;
                            tidx = test_idx == ibin;
                        end
                        
                        % run classifier
                        DEC = decoder_func(data(cell_idx(cell_num(icell,:)),idx)',groups(idx)');
                        [pre, sc] = predict(DEC,test_data(cell_idx(cell_num(icell,:)),tidx)'); % test decoder
                        
                        % Assign performance data into bins
                        p = pre;
                        r = pre;
                        for igroup = 1:group_num
                            if group_num<3;group_idx = 1;else group_idx=igroup;end
                            P{igroup}(icell,test_data_idx(tidx & test_groups==igroup)) = p(test_groups(tidx)==igroup);
                            R{igroup}(icell,data_shfl_idx(tidx & test_shfl_groups==igroup)) = ...
                                r(test_shfl_groups(tidx)==igroup);
                            S{igroup}(icell,test_data_idx(tidx & test_groups==igroup)) = sc(test_groups(tidx)==igroup,1);
                            D{igroup}{icell,ibin}(1,:) = DEC.BinaryLearners{group_idx}.Beta;
                            D{igroup}{icell,ibin}(2,:) = repmat(DEC.BinaryLearners{group_idx}.Bias,sum(cell_num(icell,:)),1);
                            D{igroup}{icell,ibin}(3,:) = mdata(cell_idx(cell_num(icell,:)));
                            D{igroup}{icell,ibin}(4,:) = sdata(cell_idx(cell_num(icell,:)));
                        end
                    end
                    Cells{irep}{icell} = cell_idx(cell_num(icell,:));
                end
                PP{irep} = P;
                RR{irep} = R;
                SC{irep} = S;
                DC{irep} = D;
            end
            
            % convert {reps}{obj}[cells trials] to {obj}[reps trials cells]
            PP = cellfun(@cell2mat,mat2cell(cellfun(@(x) permute(x,[3 2 1]),permute(reshape([PP{:}],...
                length(Data),repetitions),[2 1]),'uni',0),repetitions,ones(1,length(Data))),'uni',0);
            RR = cellfun(@cell2mat,mat2cell(cellfun(@(x) permute(x,[3 2 1]),permute(reshape([RR{:}],...
                length(Data),repetitions),[2 1]),'uni',0),repetitions,ones(1,length(Data))),'uni',0);
            SC = cellfun(@cell2mat,mat2cell(cellfun(@(x) permute(x,[3 2 1]),permute(reshape([SC{:}],...
                length(Data),repetitions),[2 1]),'uni',0),repetitions,ones(1,length(Data))),'uni',0);
            DC = num2cell(permute(reshape([DC{:}],length(Data),repetitions),[2 1]),2);
            
        end

        function [Data, Stims, info, Unit_ids] = getData(self,key,bin,stim_split)
            
            if nargin<3 || isempty(bin)
                bin = fetch1(obj.DecodeOpt & key, 'binsize');
            end
            
            % get traces
            [Traces, caTimes, keys] = getAdjustedSpikes(fuse.ActivityTrace &...
                (anatomy.AreaMembership & key),'soma'); % [time cells]
            Unit_ids = [keys.unit_id];
            
            % get rid of nans
            notnanidx = ~isnan(mean(Traces,2)); % faster than all
            Traces = Traces(notnanidx,:);
            caTimes = caTimes(notnanidx);
            if bin<0
                caTimes = caTimes - bin/1000;
                bin = abs(bin);
            end
            
            % interpolate over time
            X = @(t) interp1(caTimes-caTimes(1), Traces, t, 'linear', 'extrap');  % traces indexed by time
            
            % fetch stimuli without repeats
            trial_obj = stimulus.Trial &  ...
                ((stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"')) - ...
                (aggr(stimulus.Clip , stimulus.Trial & key, 'count(*)->n') & 'n>1')) & key;
            [flip_times, trial_idxs] = fetchn(...
                trial_obj,'flip_times','trial_idx','ORDER BY trial_idx');
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            trial_idxs = trial_idxs(tidx);
            flip_times = cell2mat(flip_times(tidx));
            Stims = unique(fetchn(stimulus.Clip & trial_obj,'movie_name'));
            
            % subsample traces
            Traces = permute(X(flip_times - caTimes(1)),[2 3 1]);
            if abs(bin)>0
                fps = 1/median(diff(flip_times(1,:)));
                d = max(1,round(abs(bin)/1000*fps));
                Traces = convn(Traces,ones(d,1)/d,'same');
                Traces = Traces(1:d:end,:,:);
            end

            Traces = permute(Traces,[2 1 3]); % in [cells bins trials]
            
            if nargin>3 && stim_split
                Data = Traces;
                info = [];
            else
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
        end
        
        function [train_groups, test_groups] = getGroups(obj,Stims,train_set,test_set)
            
            if isempty(train_set) % take all pairwise combinations of stimuli
                Stims = combnk(Stims,2);
                for igr = 1:size(Stims,1)
                    for istim = 1:size(Stims,2)
                        train_groups{igr}{istim} = Stims(igr,istim);
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
               
        function plotMasks(self,varargin)
            
            params.fontsize = 10;
            params.mi = 1;
            params.target_cell_num = [];
            params.mx = [];
            params.mn = [];
            params.tinybar = true;
            params.colormap = [];
            params.sitenum = true;
            
            params = getParams(params,varargin);
            
            % adjust colors
            if isempty(params.colormap)
%                 params.colormap = parula(30);
                colors = cbrewer('seq','YlPuBl',150);
                params.colormap = colors(1:end-40,:);
            end
            
            % get data
            [area, keys] = fetchn(self & 'brain_area <> "unknown"',...
                'brain_area');
            
            areas = unique(area);
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                idx = (strcmp(area,areas(iarea)));
                p = getPerformance(self & keys(idx),params);
                if iscell(p);p = cellfun(@nanmean,p);end
                MI{iarea} = p;
            end
            
            % plot
            figure;
            plotMask(anatomy.Masks)
            colormap(params.colormap)
            c = colorbar;
            
            % set min/max
            mx = max(cellfun(@nanmean,MI));
            mn = min(cellfun(@nanmean,MI));
            if ~isempty(params.mx);mx = params.mx;end
            if ~isempty(params.mn);mn = params.mn;end
            
            for iarea = 1:length(areas)
                mi = nanmean(MI{iarea});
                if isnan(mi);continue;end
                idx = double(uint8(floor(((mi-mn)/(mx - mn))*0.99*size(params.colormap,1)))+1);
                if params.sitenum; n = sum(~isnan(MI{iarea}));else n = [];end
                plotMask(anatomy.Masks & ['brain_area="' areas{iarea} '"'],params.colormap(idx,:),n)
            end
            
            if params.tinybar; ml = 2;else ml =5;end
            if params.mi
                labl = 'Mutual Information (bits)';
                set(c,'ytick',linspace(0,1,ml),'yticklabel',roundall(linspace(mn,mx,ml),10^-round(abs(log10(mn/10)))))
            else
                labl = 'Classification performance (%)';
                set(c,'ytick',linspace(0,1,ml),'yticklabel',round(linspace(mn*100,mx*100,ml)))
            end
            if params.tinybar; set(c,'Position',[0.88 0.14 0.025 0.3]);end
            ylabel(c,labl,'Rotation',-90,'VerticalAlignment','baseline','fontsize',params.fontsize)
            set(gcf,'name','Decoding performance - All areas')
        end
        
        function paramsExp = plotCells(self, varargin)
            params.mx_cell = [];
            params.mi = true;
            params.figure = [];
            params.fontsize = 10;
            params.linewidth = 2;
            params.colors  = [];
            params.areas = [];
            params.errors = true;
            
            params = getParams(params,varargin);
            
            [areas, trials_info, keys] = fetchn(self,'brain_area','trial_info');
            MI = getPerformance(self,params);
            cell_num = cellfun(@(x) cellfun(@length,x),cellfun(@(x) x.units{1,end},trials_info,'uni',0),'uni',0);
            
            if nargin>1 && ~isempty(params.mx_cell)
                idx = cellfun(@(x) any(x>=params.mx_cell),cell_num);
                if ~idx; return;end
                areas = areas(idx);
                MI = MI(idx);
                cell_num = cell_num(idx);
                cell_idx = repmat(find(cell_num{1}>=params.mx_cell,1,'first'),length(MI),1);
            else
                cell_idx = cellfun(@(x) length(x), cell_num);
            end
            
            if isempty(params.figure)
                figure
            end
            
            if isempty(params.areas)
                un_areas = unique(areas);
            else
                un_areas = params.areas;
            end
            
            if isempty(params.colors)
                params.colors = hsv(length(un_areas));
            end
            
            h = [];
            for iarea = 1:length(un_areas)
                area_idx = find(strcmp(areas,un_areas{iarea}));
                if isempty(params.mx_cell) 
                    for idx = area_idx(:)'
                        mi = MI(idx);
                        h(iarea) = plot([cell_num{idx}(1:cell_idx(idx))],...
                            cell2mat(cellfun(@(x) double(x(:,1:cell_idx(idx))),mi,'uni',0)'),...
                            'color',params.colors(iarea,:),'linewidth',params.linewidth);
                        hold on
                    end
                elseif ~params.errors
                    mi = MI(strcmp(areas,un_areas{iarea}));
                    h(iarea) = plot([cell_num{1}(1:cell_idx(iarea))],...
                        nanmean(cell2mat(cellfun(@(x) double(x(:,1:cell_idx(iarea))),mi,'uni',0))),...
                        'color',params.colors(iarea,:),'linewidth',params.linewidth);
                    hold on
                else
                    mi = MI(strcmp(areas,un_areas{iarea}));
                    h(iarea) = errorPlot([cell_num{1}(1:cell_idx(iarea))],...
                        cell2mat(cellfun(@(x) double(x(:,1:cell_idx(iarea))),mi,'uni',0)),...
                        'errorColor',params.colors(iarea,:),'linewidth',params.linewidth);
                end
            end
            params.l = legend(h,un_areas');
            if ~params.mi
                ylabel('Performance (% correct)')
                plot(get(gca,'xlim'),[0.5 0.5],'-.','color',[0.5 0.5 0.5]);
            elseif params.mi
                ylabel('Mutual Information (bits)')
            end
            xlabel('# of Cells')
            set(params.l,'box','off','location','northwest','fontsize',params.fontsize)
            set(gca,'fontsize',params.fontsize)
            
            if nargout; paramsExp = params;end
        end
        
        function [perf, keys] = getPerformance(self,varargin)
            
            params.target_cell_num = [];
            params.perf = 'p';
            params.mi = 1;
            params.autoconvert = true;
            
            params = getParams(params,varargin);
            
            [p, ti, keys] = fetchn(self, params.perf,'trial_info');
            
            % remove empty classes
            p = cellfun(@(x) x(~cellfun(@isempty,x)),p,'uni',0);
            perf = cell(length(p),1);
            for ikey = 1:length(p)
                if isempty(params.target_cell_num)
                    ncel =1:size(p{ikey}{1},3);
                elseif params.target_cell_num==0
                    ci = cellfun(@(x) max([find(cell2mat(cellfun(@length,x.units{1},'uni',0))>...
                        params.target_cell_num,1,'last'),0]),ti);
                    ncel = ci(ikey);
                    if ncel==0;continue;end
                else
                    ci = cellfun(@(x) max([find(cell2mat(cellfun(@length,x.units{1},'uni',0))>...
                        params.target_cell_num-1,1,'first'),0]),ti);
                    ncel = ci(ikey);
                    if ncel==0; perf{ikey} = nan;continue;end
                end
                for icell = ncel
                    P = cellfun(@(x) x(:,:,icell),p{ikey},'uni',0);
                    CM = nan(length(P));
                    for igroup = 1:length(P)
                        for itest = 1:length(P)
                            CM(igroup,itest) = nansum(P{igroup}(:)==itest);
                        end
                    end
                    if params.mi
                        perf{ikey}(end+1) = self.getMI(CM);
                    else
                        perf{ikey}(end+1) = sum(diag(CM))/sum(CM(:));
                    end
                end
            end
            if all(cellfun(@length,perf)==1) && params.autoconvert
                perf = cell2mat(perf);
            end
        end
        
        function params = getParams(self, param_type)
            [trial_info, binsize] = fetch1(self*obj.DecodeOpt,'trial_info','binsize');
            for iclass = 1:length(trial_info.clips)
                keys = cell2struct([num2cell(trial_info.clips{iclass},2),...
                    trial_info.names{iclass}']',{'clip_number','movie_name'},1);
                params{iclass} = getParam(stimulus.MovieParams, keys, param_type, ...
                    trial_info.bins{iclass}*abs(binsize)/1000 - abs(binsize)/2/1000);
            end
        end
        
        function cell_num = getCellNum(self, key)
            if nargin<2; restr = proj(self);else restr = key;end
            ti = fetchn(obj.Dec & restr,'trial_info');
            cell_num = [];
            for ikey = 1:length(ti)
                cell_num{ikey} = cellfun(@length,ti{ikey}.units{1});
            end
            if ikey==1;cell_num=cell_num{1};end
        end
        
        function plot(self,varargin)
            
            params.fontsize = 10;
            params.mi = 1;
            params.target_cell_num = [];
            params.mx = [];
            params.mn = [];
            params.tinybar = true;
            params.colormap = [];
            params.type = 'boxfun';
            
            params = getParams(params,varargin);
            
            % adjust colors
            if isempty(params.colormap)
                params.colormap = parula(30);
            end
            
            % get data
            [area, keys] = fetchn(self & 'brain_area <> "unknown"',...
                'brain_area');
            
            areas = unique(area);
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                idx = (strcmp(area,areas(iarea)));
                p = getPerformance(self & keys(idx),params);
                if iscell(p);p = cellfun(@nanmean,p);end
                MI{iarea} = p;
            end
            
            switch params.type
                case 'boxfun'
                    boxfun(MI,'barwidth',0.9,'names',areas,'sig',1,'rawback',1)
                case 'barfun'
                    barfun(MI,'barwidth',0.9,'names',areas,'sig',1)
                case 'masks'
                    plotMasks(self,params)
            end
        end
    end
    
    methods (Static)
        function mi = getMI(CM)
            CM = CM+eps;
            p = CM/sum(CM(:));
            pi = sum(CM,2)/sum(CM(:));
            pj = sum(CM,1)/sum(CM(:));
            pij = pi*pj;
            mi = sum(sum(p.*log2(p./pij)));
        end
    end
end
