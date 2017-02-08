%{
mov3d.Repeats (computed) # calcium trace
-> preprocess.Sync
-> mov3d.RepeatsOpt
-> preprocess.SpikeMethod
-> preprocess.Method
---
r                    : longblob                      # reliability
%}

classdef Repeats < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  = (experiment.Scan  ...
            * (preprocess.Spikes))...
            * (mov3d.RepeatsOpt & 'process = "yes"') ...
            * (preprocess.Sync  & (vis.MovieClipCond * vis.Trial & ...
            (vis.Movie & 'movie_class="object3d" OR movie_class = "multiobjects"') & ...
            'trial_idx between first_trial and last_trial'))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            method = fetch1(mov3d.RepeatsOpt & key,'method');
            
            Data = getData(self,key); % [Cells, Time, Trials]
            
            switch method
                case 'explainedVar'
                    r = mean(cellfun(@reliability,Data));
                case 'corr'
                    r = nan(length(Data),1);
                    for iobj = 1:length(Data)
                        cor = nan(size(Data{iobj},1),1);
                        for icell = 1:size(Data{iobj},1)
                           traceZ = zscore(squeeze(Data{iobj}(icell,:,:)));
                           c = corr(traceZ);
                           cor(icell) = nanmean(c(logical(tril(ones(size(traceZ,2)),-1))));
                        end
                        r(iobj) = nanmean(cor);
                    end
                    r = nanmean(r);
            end
            % insert
            tuple.r = r; 
            self.insert(tuple)
            
        end
    end
    
    methods
        function Data = getData(obj,key,ibin) % {uni_stims}(cells,time,repeats)
            
            [bin] = fetch1(mov3d.RepeatsOpt & key, 'binsize');
            if nargin>2;bin = ibin;end
            
            % get stuff
            [Traces, caTimes] = pipetools.getAdjustedSpikes(key);
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            [flip_times, trial_idxs, mov, clip] = (fetchn(trials * vis.MovieClipCond,'flip_times','trial_idx','movie_name','clip_number'));
            
            % filter out incomplete trials
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
  
            % find repeated trials
            unimov=unique(mov);
            stim_index = repmat((1:length(unimov))',1,size(mov,1));
            stim_index = stim_index(strcmp(repmat(unimov,1,size(mov,1)),repmat(mov',length(unimov),1)));
            stimuli = [stim_index clip];
            [~,~,uni_idx] = unique(stimuli,'rows');
            ridx = arrayfun(@(x) sum(x==uni_idx),uni_idx)>1;
            unis = unique(stimuli(ridx,:),'rows');
            
            % Subsample traces
            flip_times = cell2mat(flip_times(tidx & ridx));
            xm = min([length(caTimes) length(Traces)]);
            X = @(t) interp1(caTimes(1:xm)-caTimes(1), Traces(1:xm,:), t, 'linear', nan);  % traces indexed by time
            fps = 1/median(diff(flip_times(1,:)));
            d = max(1,round(bin/1000*fps));
            traces = convn(permute(X(flip_times - caTimes(1)),[2 3 1]),ones(d,1)/d,'same');
            traces = traces(1:d:end,:,:);
            
            % split for unique stimuli
            Data = [];k=[];
            for istim = 1:length(unis)
                k.movie_name = unimov{unis(istim,1)};
                k.clip_number = unis(istim,2);
                stim_trials = fetchn(trials*vis.MovieClipCond & k,'trial_idx');
                Data{istim} = permute(traces(:,:,ismember(trial_idxs(tidx & ridx),stim_trials)),[2 1 3]);
            end
        end
    end
end
