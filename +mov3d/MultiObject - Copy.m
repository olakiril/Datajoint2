%{
mov3d.MultiObject (computed) # calcium trace
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

classdef MultiObject < dj.Relvar & dj.AutoPopulate
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
            
            angle = @(x) acosd(sum(x)./sqrt(sum(x.^2))/sqrt(size(x,1))); % x in [cells trials]
            
            [decoder,k_fold,shuffle,train_set,test_set,repetitions] = ...
                fetch1(mov3d.DecodeMultiOpt & key,...
                'decoder','k_fold','shuffle','train_set','test_set','repetitions');
            
            [Traces, Stims] = getData(self,key); % [Cells, Obj, Trials]
     
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
            
            bin = fetch1(mov3d.DecodeMultiOpt & key, 'binsize');

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

            % get Traces for all common clips
            Data = [];
            [multi_stim, stim_cla] = getGroups(obj,Stims);
            uni_stims = unique(stim_cla);
            for istim = 1:length(uni_stims)
                stim_trials = [];clip_num = [];
                singlestims = Stims(strcmp(stim_cla,uni_stims{istim}) & ~ multi_stim);
                for iobj = 1:length(singlestims)
                    key.movie_name = singlestims{iobj};
                    [stim_trials{iobj}, clip_num{iobj}] = fetchn(trials*vis.MovieClipCond & key,'trial_idx','clip_number');
                end
                [uni_clips,uni_idx{1},uni_idx{2}] = intersect(clip_num{1},clip_num{2});
                for iobj = 1:length(singlestims)
                    Data{strcmp(singlestims{iobj},Stims)} = reshape(permute(traces(:,:,ismember(trial_idxs,stim_trials{iobj}(uni_idx{iobj}))),[2 1 3]),size(traces,2),[]);
                end

                key.movie_name = Stims{strcmp(stim_cla,uni_stims{istim}) & multi_stim};
                [stim_trials, clips] = fetchn(trials*vis.MovieClipCond & key,'trial_idx','clip_number');
                stim_trials = stim_trials(ismember(clips,uni_clips));
                Data{strcmp(key.movie_name,Stims)} = reshape(permute(traces(:,:,ismember(trial_idxs,stim_trials)),[2 1 3]),size(traces,2),[]);
            end
        end
        
        function [multi_stim, stim_class] = getGroups(obj,Stims)
            multi_stim = cell2mat(cellfun(@length,regexp(Stims,'(?<=(o))\d+','match'),'uni',0))>1;
            stim_class = cellfun(@cell2mat,regexp(Stims,'(?<=(s))\d+','match'),'uni',0);
        end
        
    end
end
