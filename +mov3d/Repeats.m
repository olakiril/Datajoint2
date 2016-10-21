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
            * (preprocess.Spikes & 'spike_method = 3'  & 'extract_method=2'))...
            * (mov3d.RepeatsOpt & 'process = "yes"') ...
            * (preprocess.Sync & (vis.MovieClipCond & (vis.Movie & 'movie_class="object3d"')))
        
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
        function Data = getData(obj,key,ibin)
            
            [bin, rf_idx] = fetch1(mov3d.RepeatsOpt & key, 'binsize','restrict_rf');
            if nargin>2;bin = ibin;end
            
            if rf_idx > 0
                index = true;
                [rf_idx, rf_trials] = fetch1(mov3d.RFFilter & key,'rf_idx','rf_trials');
            else
                index = false;
            end
            
            [Traces, caTimes] = pipetools.getAdjustedSpikes(key);
            xm = min([length(caTimes) length(Traces)]);
            X = @(t) interp1(caTimes(1:xm)-caTimes(1), Traces(1:xm,:), t, 'linear', nan);  % traces indexed by time
            
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            [mov,clip] = fetchn(trials*vis.MovieClipCond,'movie_name','clip_number');
            movi = strcmp(mov,'obj2v3')+1;
            [~,~,y] = unique([movi clip],'rows');
            idx = arrayfun(@(x) sum(x==y),y)>1;
            trials = fetch(trials*vis.MovieClipCond, '*', 'ORDER BY trial_idx'); %fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx') 2016-08
            trials = trials(idx);
            unis = unique(y(idx));
            
            snippet = []; % traces: {stimulus}(subbin,cells)
            for itrial = 1:length(trials)
                trial = trials(itrial);
                
                % extract relevant trace & bin
                fps = 1/median(diff(trial.flip_times));
                t = trial.flip_times - caTimes(1);
                d = max(1,round(bin/1000*fps));
                trace = convn(X(t),ones(d,1)/d,'same');
                trace = trace(1:d:end,:);
                if index; trace = trace(rf_idx{trial.trial_idx==rf_trials},:);end
                snippet{itrial} = trace;
            end
            
            Data = [];
            for iuni = 1:length(unis)
                dat = snippet(unis(iuni)==y(idx));
                dat = dat(~cellfun(@isempty,dat));
                Data{iuni} = permute(cell2mat(cellfun(@(x) permute(x,[1 3 2]),dat,'uni',0)),[3 1 2]);
            end
        end
    end
end
