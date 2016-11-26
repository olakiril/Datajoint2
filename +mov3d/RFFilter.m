%{
mov3d.RFFilter (computed) # population RF
-> preprocess.Sync
-> preprocess.SpikeMethod
-> preprocess.Method
binsize                : float           # bin size 
---
map                    : longblob        # standard deviation from center of population RF
rf_idx                 : longblob        # index of subbins with stimulus in pop RF (trials)(subbins)
rf_trials              : longblob        # trial index (trial indexes)
v1_rfs                 : longblob        # Individual RF centers [x y]
%}

classdef RFFilter < dj.Relvar & dj.AutoPopulate
    
    properties        
         popRel  = (experiment.Scan  ...
            * (pro(preprocess.Spikes,'spike_method->spike_inference','extract_method->segment_method') ...
            & 'spike_inference = 2'  & 'segment_method=2'))...
            * (pro(mov3d.DecodeOpt & 'process = "yes"' & 'restrict_rf>0','binsize->bin')) ...
            * (preprocess.Sync & (vis.MovieClipCond & (vis.Movie & 'movie_class="object3d"'))) ...
            & tuning.MonetFit 
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = rmfield(key,'dec_opt');
            
            % params
            [binsize, rf_thr] = fetch1(mov3d.DecodeOpt & key,'binsize','restrict_rf');
            tuple.binsize = binsize;
            
            if exists(self & tuple); return;end
            
            % get V1 RFs
            sesskey = rmfield(key,'scan_idx');
            V1key = fetch(experiment.Scan & 'brain_area= "V1"' & sesskey);
            if isempty(V1key)
                sesskey.session = sesskey.session-1;
                V1key = fetch(experiment.Scan & 'brain_area= "V1"' & sesskey);
            end
            [xloc, yloc] = fetchn(tuning.MonetFit & V1key(1),'x','y');
            
            if isempty(xloc); warning('No RFs found!');return;end
            
            % convert to pixel space
            sess = fetch(preprocess.Sync*vis.Session & key,'resolution_x','resolution_y','monitor_distance','monitor_size');
            rect = [sess.resolution_x sess.resolution_y];
            degPerPix = 180/pi*sess.monitor_size*2.54/norm(rect(1:2))/sess.monitor_distance;
            xloc = xloc/degPerPix;
            yloc = yloc/degPerPix;
            
            % fit RFs
            m = fitgmdist([xloc,yloc],1);
            
            % create rf boundary
            [x,y] = meshgrid((1:rect(1))-rect(1)/2,(1:rect(2))-rect(2)/2);
            X=[x(:) y(:)];
            X = bsxfun(@minus, X, m.mu);
            d = sum((X /m.Sigma) .* X, 2);
            pop_rf = reshape(d,rect(2),rect(1)); % in pixel space
            
            % stimulus_trial_xy_position
            [paramsObj,obj,fps] = fetchn(vis.Movie & (vis.MovieClipCond & key),...
                'params','movie_name','frame_rate');
            
            stim_idx = [];
            for iobj = 1:length(obj)
                params = paramsObj{iobj};
                frameStep = fps(iobj)*binsize/1000; % in frames
                frameIdx = 1:frameStep:params.frames(end);
                px = interpn(params.frames,params.camera_pos_x,frameIdx,'cubic');
                pz = interpn(params.frames,params.camera_pos_z,frameIdx,'cubic');
                nx = round(normalize(px)*rect(1));
                nz = round(normalize(pz)*rect(2));
                nx(nx==0) = 1;
                nz(nz==0) = 1;
                
                % stimulus_bin_xy_position inside circle
                stim_idx{iobj} = pop_rf(sub2ind(size(pop_rf),nz,nx))<rf_thr;
            end
            
            % get trials
            trials = pro(preprocess.Sync*vis.Trial & (experiment.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
            trials = fetch(trials*vis.MovieClipCond, '*', 'ORDER BY trial_idx'); %fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx') 2016-08
            
            % find bins within the pop RF
            rf_idx = []; rf_trials = [];
            for itrial = 1:length(trials)
                rf_trials(itrial) = trials(itrial).trial_idx;
                obj_idx = strcmp(obj,trials(itrial).movie_name);
                frames_per_trial = trials(itrial).cut_after*fps(obj_idx);
                start = (trials(itrial).clip_number - 1)*frames_per_trial;
                rf_idx{itrial} = stim_idx{obj_idx}(find(frameIdx>start,1,'first') : ...
                    find(frameIdx<(start+frames_per_trial),1,'last'));
            end
            
            % populate
            tuple.map = pop_rf;
            tuple.rf_idx = rf_idx;
            tuple.rf_trials = rf_trials;
            tuple.v1_rfs = [xloc yloc];
            self.insert(tuple)
            
        end
    end
    
end
