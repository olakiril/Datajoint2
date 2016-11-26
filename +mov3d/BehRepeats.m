%{
mov3d.BehRepeats (computed) # calcium trace
-> preprocess.Sync
-> mov3d.BehRepeatsOpt
-> preprocess.SpikeMethod
-> preprocess.Method
---
r_moving                   : longblob                      # reliability
r_still                  : longblob                        # reliability
moving                   : longblob
%}

classdef BehRepeats < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  = (experiment.Scan  ...
            * (preprocess.Spikes & 'spike_method = 5'  & 'extract_method=2'))...
            * (mov3d.BehRepeatsOpt & 'process = "yes"') ...
            * (preprocess.Sync & (vis.MovieClipCond & (vis.Movie & 'movie_class="object3d"')))
        
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            method = fetch1(mov3d.RepeatsOpt & key,'method');
            
            [mData, sData, moving] = getData(self,key); % [Cells, Time, Trials]
            
            switch method
                case 'explainedVar'
                    r_moving = nanmean(cellfun(@reliability,mData));
                    r_still = nanmean(cellfun(@reliability,sData));
                case 'corr'
                    r = nan(length(mData),1);
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
            tuple.r_moving = r_moving; 
            tuple.r_still = r_still; 
            tuple.moving = moving;
            self.insert(tuple)
            
        end
    end
    
    methods
        function [mData, sData, moving] = getData(obj,key,ibin)
            
            [bin, rf_idx, beh_thr] = fetch1(mov3d.BehRepeatsOpt & key, 'binsize','restrict_rf','beh_thr');
            if nargin>2;bin = ibin;end
            
            if rf_idx > 0
                index = true;
                [rf_idx, rf_trials] = fetch1(mov3d.RFFilter & key,'rf_idx','rf_trials');
            else
                index = false;
            end
            
            nslices = fetch1(preprocess.PrepareGalvo & key, 'nslices');
            BehTimes = fetch1(preprocess.BehaviorSync &  (experiment.Scan & key), 'frame_times');
            BehTimes = BehTimes(1:nslices:end);

            [vel, tm] = fetch1(preprocess.Treadmill & key,'treadmill_vel','treadmill_time');
            vel = abs(vel);
            mvel = mean(vel);
            svel = std(vel);
            
            [Traces, caTimes] = pipetools.getAdjustedSpikes(key);
            xm = min([length(caTimes) length(Traces)]);
            zm = min([length(vel) length(tm)]);
            X = @(t) interp1(caTimes(1:xm)-caTimes(1), Traces(1:xm,:), t, 'linear', nan);  % traces indexed by time
            Y = @(t) interp1(caTimes(1:xm)-caTimes(1), BehTimes(1:xm,:), t, 'linear', nan);  % traces indexed by Beh time
            Z = @(t) interp1(tm(1:zm), vel(1:zm,:), t, 'linear', nan);  % vel indexed by Beh time
           
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            [mov,clip] = fetchn(trials*vis.MovieClipCond,'movie_name','clip_number');
            movi = strcmp(mov,'obj2v3')+1;
            [~,~,y] = unique([movi clip],'rows');
            idx = arrayfun(@(x) sum(x==y),y)>1;
            trials = fetch(trials*vis.MovieClipCond, '*', 'ORDER BY trial_idx'); %fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx') 2016-08
            trials = trials(idx);
            unis = unique(y(idx));
            
            snippet = []; % traces: {stimulus}(subbin,cells)
            moving = [];
            for itrial = 1:length(trials)
                trial = trials(itrial);
                
                % extract relevant trace & bin
                fps = 1/median(diff(trial.flip_times));
                t = trial.flip_times - caTimes(1);
                
                ball = Z(Y(t));
                moving(itrial) = mean(ball)>(mvel+beh_thr*svel);
                d = max(1,round(bin/1000*fps));
                trace = convn(X(t),ones(d,1)/d,'same');
                trace = trace(1:d:end,:);
                if index; trace = trace(rf_idx{trial.trial_idx==rf_trials},:);end
                snippet{itrial} = trace;
            end
            
            mData = [];
            sData = [];
            for iuni = 1:length(unis)
              
                dat = snippet(unis(iuni)==y(idx) & moving');
                dat = dat(~cellfun(@isempty,dat));
                mData{iuni} = permute(cell2mat(cellfun(@(x) permute(x,[1 3 2]),dat,'uni',0)),[3 1 2]);
                
                dat = snippet(unis(iuni)==y(idx) & ~moving');
                dat = dat(~cellfun(@isempty,dat));
                sData{iuni} = permute(cell2mat(cellfun(@(x) permute(x,[1 3 2]),dat,'uni',0)),[3 1 2]);
                
            end
        end
    end
end
