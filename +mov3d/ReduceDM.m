%{
mov3d.ReduceDM (computed) # calcium trace
-> mov3d.ReduceDMOpt
-> preprocess.Sync
-> preprocess.SpikeMethod
-> preprocess.Method
---
mapped                      : longblob             # reduced data
params                      : longblob             # stimulus parameters
trials                      : longblob             # trial idx
labels                      : longblob             # object idx
%}

classdef ReduceDM < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  = (experiment.Scan  ...
            * (preprocess.Spikes & 'spike_method = 3'  & 'extract_method=2'))...
            * (mov3d.ReduceDMOpt & 'process = "yes"') ...
            * (preprocess.Sync & (vis.MovieClipCond & (vis.Movie & 'movie_class="object3d"')))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            % setup params
            [bin,no_dims,initial_dims,perplexity,gwin] = ...
                fetch1(mov3d.ReduceDMOpt & key,...
                'binsize','dimensions','initial_dims','perplexity','gauss_win');
            gwin = floor(gwin/bin); % convert to time bins
            
            % get data
            [Data, Trials, Params] = getData(self,key,bin);
            
            % get rid of Data
            [~,~,idx] = ind2sub(size(Data),find(isnan(Data)));
            Data(:,:,idx) =[];
            Trials(:,idx) = [];
            Params(:,:,idx) = [];
            
            % arrange data
            Labels = nan(size(Trials));
            Labels(1,:) = 1;
            Labels(2,:) = 2;
            Labels = Labels(:);
            trials = Trials(:);
            uTrials = unique(trials);
            NewData = Data(:,:);
            NewParams = Params(:,:);
            
            % prefilter data
            for itrial = 1:length(uTrials);
                trialIdx = trials==uTrials(itrial);
                for icell = 1:size(Data,1)
                    NewData(icell,trialIdx) = conv(NewData(icell,trialIdx),gausswin(gwin),'same');
                end
            end
            
            % do it
            mappedX = tsne(NewData', [], no_dims, initial_dims, perplexity);
            
            % correct for key mismach
            tuple = rmfield(tuple,'spike_method');
            tuple = rmfield(tuple,'extract_method');
            tuple.spike_inference = key.spike_method;
            tuple.segment_method = key.extract_method;
            
            % insert
            tuple.mapped = mappedX;
            tuple.params = NewParams;
            tuple.trials = trials;
            tuple.labels = Labels;
            self.insert(tuple)
            
        end
    end
    
    methods
        function [Data, Trials, Params] = getData(obj,key,ibin)
            
            if nargin<3
                bin = fetch1(mov3d.ReduceDMOpt & key, 'binsize');
            else
                 bin = ibin;
            end
            
            [Traces, caTimes] = pipetools.getAdjustedSpikes(key);
            xm = min([length(caTimes) size(Traces,1)]);
            X = @(t) interp1(caTimes(1:xm)-caTimes(1), Traces(1:xm,:), t, 'linear', nan);  % traces indexed by time
            
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            trials = fetch(trials*vis.MovieClipCond, '*', 'ORDER BY trial_idx'); %fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx') 2016-08
            
            snippet = []; % traces: {object,trials}(subbin,cells)
            stims = [2 1];
            idx = 0;
            trial_idx = [];
            for trial = trials'
                idx = idx+1;
                stim = stims(~isempty(strfind(trial.movie_name,'obj1'))+1); % stims(~isempty(strfind(trial.path_template,'obj1'))+1); 2016-08
                % extract relevant trace & bin
                fps = 1/median(diff(trial.flip_times));
                t = trial.flip_times - caTimes(1);
                d = max(1,round(bin/1000*fps));
                trace = convn(X(t),ones(d,1)/d,'same');
                trace = trace(1:d:end,:);
                snippet{stim,idx} = trace;
                trial_idx{stim,idx} = repmat(trial.trial_idx,size(trace,1),1);
            end
            
            A = snippet(1,:);
            A_trials = trial_idx(1,:);
            idx = ~cellfun(@isempty,A);
            A = A(idx);
            A_trials = A_trials(idx);
            objA = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),A,'uni',0)'),size(A{1},2),[]),[3 1 2]);
            objA_trials = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),A_trials,'uni',0)'),size(A_trials{1},2),[]),[3 1 2]);
            
            B = snippet(2,:);
            B_trials = trial_idx(2,:);
            idx = ~cellfun(@isempty,B);
            B = B(idx);
            B_trials = B_trials(idx);
            objB = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),B,'uni',0)'),size(B{1},2),[]),[3 1 2]);
            objB_trials = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),B_trials,'uni',0)'),size(B_trials{1},2),[]),[3 1 2]);
            
            % get params
            [params, param_trials] = getParams(obj,key,bin);
            objA_params = cell2mat(params(ismember(param_trials,cellfun(@(x) x(1),A_trials)))');
            objB_params = cell2mat(params(ismember(param_trials,cellfun(@(x) x(1),B_trials)))');
            
            % Arrange data
            mS = min([size(objA,3) size(objB,3)]);
            Data = reshape(permute(objA(:,:,1:mS),[2 4 3 1]),size(objA,2),1,[]);
            Data(:,2,:) = reshape(permute(objB(:,:,1:mS),[2 4 3 1]),size(objB,2),1,[]);
            Trials = reshape(permute(objA_trials(:,:,1:mS),[2 4 3 1]),size(objA_trials,2),1,[]);
            Trials(:,2,:) = reshape(permute(objB_trials(:,:,1:mS),[2 4 3 1]),size(objB_trials,2),1,[]);
            Trials = squeeze(Trials(1,:,:));
            Params = permute(objA_params(1:mS,:),[2 3 1]);
            Params(:,2,:) = permute(objB_params(1:mS,:),[2 3 1]);
            
        end
        
        function [params, param_trials] = getParams(obj,key,bin)
            
            speed = @(x,y,timestep) sqrt(x.^2+y.^2)./timestep;
            binsize= fetch1(mov3d.ReduceDMOpt & key, 'binsize');
            if nargin>2;binsize = bin;end
            
            % stimulus_trial_xy_position
            [paramsObj,obj,fps] = fetchn(vis.Movie & (vis.MovieClipCond & key),...
                'params','movie_name','frame_rate');
            
            int_params = [];
            for iobj = 1:length(obj)
                timestep = mean(diff(paramsObj{iobj}.frames))/fps(iobj);
                par = struct2array(paramsObj{iobj});
                
                par = par(:,[14 16:end]);
                par(:,end+1) = [0;speed(diff(par(:,1)),diff(par(:,2)),timestep)];
                frameStep = fps(iobj)*binsize/1000; % in frames
                frameIdx = 1:frameStep:paramsObj{iobj}.frames(end);
                int_params{iobj} = nan(length(frameIdx),size(par,2));
                for iparam = 1:size(par,2)
                    int_params{iobj}(:,iparam) = interpn(paramsObj{iobj}.frames,par(:,iparam),frameIdx,'cubic');
                end
            end
            
            % get trials
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            trials = fetch(trials*vis.MovieClipCond, '*', 'ORDER BY trial_idx'); %fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx') 2016-08
            
            % find bins within the pop RF
            params = []; param_trials = [];
            
            for itrial = 1:length(trials)
                param_trials(itrial) = trials(itrial).trial_idx;
                obj_idx = strcmp(obj,trials(itrial).movie_name);
                frames_per_trial = trials(itrial).cut_after*fps(obj_idx);
                start = (trials(itrial).clip_number - 1)*frames_per_trial;
                params{itrial} = int_params{obj_idx}(find(frameIdx>start,1,'first') : ...
                    find(frameIdx<(start+frames_per_trial),1,'last'),:);
            end
        end
        
        function plotParams(obj,param,object)
            
            [NewParams,trials,mappedX,Labels] = fetch1(obj,'params','trials','mapped','labels');
            
            if nargin<2;param = 1;end
            if nargin<3;object=1;end
            
            tbin = 3;
            uTrials = unique(trials);
%             figure
%             hold on;
            param = squeeze(normalize(NewParams(param,:)));
            
            for itrial = 2:length(uTrials)
                
                trialIdx = trials==uTrials(itrial);
                tIdx = all(Labels(trialIdx)==object);
                
                if tIdx;continue;end
                
                
                xx = interpn(mappedX(trialIdx,1),tbin,'cubic');
                yy = interpn(mappedX(trialIdx,2),tbin,'cubic');
                zz = interpn(mappedX(trialIdx,3),tbin,'cubic');
                
                col = interpn(param(trialIdx)',tbin,'cubic');
                
                surface([xx; xx],[yy; yy],[zz; zz],[col;col],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',1,'facealpha',0.5,'edgealpha',0.5);
                
                
            end
            
            colormap jet
        end
    end
    
end
