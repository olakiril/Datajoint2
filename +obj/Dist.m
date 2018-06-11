%{
# Distance between objects
-> fuse.ScanDone
-> anatomy.Area
-> obj.DistOpt
---
-> stimulus.Sync
distance                 : mediumblob                      # avg distance between stimuli
stimuli                  : mediumblob                      # stimuli identity
%}

classdef Dist < dj.Computed
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.DistOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"'))))
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % get train & test groups
            [obj_set,repetitions,units] = fetch1(obj.DistOpt & key,'obj_set','repetitions','units');
            
            % get DAta
            [Traces, Stims] = getData(self,key); % [Cells, Obj, Trials]
            obj_groups = getGroups(self,Stims,obj_set);
            
            % run the decoding
            Dist = nan(length(obj_groups),repetitions);
            for iGroup = 1:length(obj_groups)
                obj_data = [];
                for iClass = 1:length(obj_groups{iGroup})
                    tgroup = obj_groups{iGroup}{iClass};
                    stim_idx = any(strcmp(...
                        repmat(tgroup,size(Stims,1),1),...
                        repmat(Stims,1,size(tgroup,2)))',1);
                    obj_data{iClass} = cell2mat(Traces(stim_idx));
                end
                
                mn_sz = min(cellfun(@(x) size(x,2),obj_data));
                for irep = 1:repetitions
                    unit_idx = randperm(size(obj_data{1},1),units);
                    stim_idx1 = randperm(size(obj_data{1},2),mn_sz);
                    stim_idx2 = randperm(size(obj_data{2},2),mn_sz);
                    dist = pdist2(obj_data{1}(unit_idx,stim_idx1)',obj_data{2}(unit_idx,stim_idx2)');
                    Dist(iGroup,irep) = nanmean(dist(:));
                end
            end
            
            % insert
            key.distance = Dist;
            key.stimuli = obj_groups;
            self.insert(key)
        end
    end
    
    methods
        function [Data, Stims, info, Unit_ids] = getData(self,key,bin,stim_split)
            
            if nargin<3 || isempty(bin)
                bin = fetch1(obj.DistOpt & key, 'binsize');
            end
            
            % get traces
            [Traces, caTimes, keys] = getAdjustedSpikes(fuse.ActivityTrace & (anatomy.AreaMembership & key),'soma'); % [time cells]
            Unit_ids = [keys.unit_id];
            
            % get rid of nans
            notnanidx = ~isnan(mean(Traces,2)); % faster than all
            Traces = Traces(notnanidx,:);
            caTimes = caTimes(notnanidx);
            
            % interpolate over time
            X = @(t) interp1(caTimes-caTimes(1), Traces, t, 'linear', 'extrap');  % traces indexed by time
            
            % fetch stimuli without repeats
            [flip_times, trial_idxs] = fetchn(...
                stimulus.Trial &  ...
                ((stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"')) - ...
                (aggr(stimulus.Clip , stimulus.Trial & key, 'count(*)->n') & 'n>1')) & key,...
                'flip_times','trial_idx','ORDER BY trial_idx');
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            trial_idxs = trial_idxs(tidx);
            flip_times = cell2mat(flip_times(tidx));
            Stims = unique(fetchn(stimulus.Clip &  (stimulus.Trial & key),'movie_name'));
            
            % subsample traces
            Traces = permute(X(flip_times - caTimes(1)),[2 3 1]);
            if bin>0
                fps = 1/median(diff(flip_times(1,:)));
                d = max(1,round(bin/1000*fps));
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
        
        function [train_groups] = getGroups(obj,Stims,train_set)
            
            if isempty(train_set) % take all pairwise combinations of stimuli
                Stims = combnk(Stims,2);
                for igroup = 1:size(Stims,1)
                    for istim = 1:size(Stims,2)
                        train_groups{igroup}{istim} = Stims(igroup,istim);
                    end
                end
            else
                train_groups = splitGroups(train_set);
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
        
                function plotMasks(self)
            
            % get data
            [perf, area] = fetchn(self & 'brain_area <> "unknown"',...
                'distance','brain_area');
            areas = unique(area);
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                idx = find(strcmp(area,areas(iarea)));
               
                    labl = 'Avg. Eucl. Distance';
                    mi = nan(length(idx),1);
                    for iscan = 1:length(idx)
                         mi(iscan) = nanmean(reshape(perf{idx(iscan)},1,[]));
                    end
                    MI{iarea} = mi;
            end
            
            % plot
            f = figure;
            colors = parula(30);
            plotMask(anatomy.Masks)
            colormap parula
            c = colorbar;
            ylabel(c,labl,'Rotation',-90,'VerticalAlignment','baseline')
            
            
            mx = max(cellfun(@nanmean,MI));
            mn = min(cellfun(@nanmean,MI));
            for iarea = 1:length(areas)
                mi = nanmean(MI{iarea});
                if isnan(mi);continue;end
                idx = double(uint8(floor(((mi-mn)/(mx - mn))*0.99*size(colors,1)))+1);
                plotMask(anatomy.Masks & ['brain_area="' areas{iarea} '"'],colors(idx,:),sum(~isnan(MI{iarea})))
            end
            if nargin>1 && norm
            set(c,'ytick',linspace(0,1,5),'yticklabel',roundall(linspace(mn,mx,5),mx/10))
            else
                set(c,'ytick',linspace(0,1,5),'yticklabel',round(linspace(mn*100,mx*100,5)))
            end
        end
        

    end
    
end
