%{
# Site statistics for object stimuli
-> fuse.ScanDone
-> anatomy.Area
-> obj.SiteStatsOpt
-> stimulus.MovieClass
---
-> stimulus.Sync
neurons               : smallint           #
distIn                : float              # average eucleadian distance
distTrans             : float              # across object identity average eucleadian distance
distCis               : float              # within object identity average eucleadian distance
corr                  : float              #
mean                  : float              # mean
variance              : float              # variance
tempwin               : float              # response temporal window
pspars                : float              # population sparseness
pzero                 : float              # probability of zero neurons responding
pkurt                 : float              # kurtosis
lspars                : float              # lifetime sparseness
lzero                 : float              # probability of zero response
lkurt                 : float              # lifetime kurtosis
%}

classdef SiteStats < dj.Computed
    
    properties
        keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.SiteStatsOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"'))))
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % get Data
            [Traces, Stims] = getData(self,key); % [Cells, Obj, Trials]
            ObjID = cellfun(@(x) str2num(x{1}),regexp(Stims,'\d(?=\w*v\d)','match'));
            
            % Loop through stimuli
            for iStim = 1:length(Stims)
                traces = Traces{iStim};
                
                key.neurons = size(traces,1);
                key.mean = nanmean(traces(:));
                key.variance = mean(nanvar(traces,[],2));
                key.tempwin = calcResponseWindow(traces');
                c =  corr(traces');
                key.corr = nanmean(c(logical(tril(ones(size(c)),-1))));
                key.distIn = nanmean(pdist(traces));
                group = ObjID==ObjID(iStim);
                key.distTrans = nanmean(reshape(pdist2(traces',cell2mat(Traces(~group))'),[],1));
                group(iStim) = false;
                key.distCis = nanmean(reshape(pdist2(traces',cell2mat(Traces(group))'),[],1));
                key.pspars = nanmean(sparseness(traces'));
                key.pzero = nanmean(sparseness(traces','type','pzero'));
                key.pkurt = nanmean(sparseness(traces','type','kurtosis'));
                key.lspars = nanmean(sparseness(traces));
                key.lzero = nanmean(sparseness(traces,'type','pzero'));
                key.lkurt = nanmean(sparseness(traces,'type','kurtosis'));
                
                % insert
                key.movie_name = Stims{iStim};
                self.insert(key)
                
            end
        end
    end
    
    methods
        function [Data, Stims, info, Unit_ids] = getData(self,key,bin,stim_split)
            
            if nargin<3 || isempty(bin)
                bin = fetch1(obj.DecodeOpt & key, 'binsize');
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
    end
end
