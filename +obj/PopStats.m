%{
-> stimulus.Sync
-> fuse.ScanDone
-> anatomy.Area
-> obj.SiteStatsOpt
-> stimulus.Movie
---
dist_in                     : float                        # average eucleadian distance
dist_trans=null             : float                        # across object identity average eucleadian distance
dist_cis=null               : float                        # within object identity average eucleadian distance
dist_mean                   : float                        # between spherical means
corr                        : float                        # average total correlation across neurons
smean                       : float                        # spherical mean - Angle from diagonal
svar                        : float                        # spherical variance
rmean                       : float                        # radial mean 
rvar                        : float                        # radial variance
%}

classdef PopStats < dj.Computed
    
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
            [neurons, reps] = fetch1(obj.SiteStatsOpt & key,'neurons','reps');
            assert(neurons < size(Traces{1},1),sprintf('Site has only %d neurons',size(Traces{1},1)))
            
            % Loop through stimuli
            for iStim = 1:length(Stims)
                
                [corrs, dist_in, dist_trans, dist_cis, dist_mean, smean, svar, rmean, rvar] = initialize('nan',reps,1);
                for irep = 1:reps
                    rseed = RandStream('mt19937ar','Seed',irep);
                    neuro_idx = randperm(rseed,size(Traces{iStim},1),neurons);
                    traces = Traces{iStim}(neuro_idx,:);
                    c =  corr(traces');
                    corrs(irep) = nanmean(c(logical(tril(ones(size(c)),-1))));
                    dist_in(irep) = nanmean(pdist(traces'));
                    group = ObjID==ObjID(iStim);
                    traces_trans = cell2mat(cellfun(@(x) x(neuro_idx,:),Traces(~group),'uni',0));
%                     dist_trans(irep) = nanmean(reshape(pdist2(traces',traces_trans(:,randperm(rseed,size(traces_trans,2),size(traces,2)))'),[],1));
                    dist_trans(irep) = mean(pdist2(traces',mean(traces_trans,2)'));

                    group(iStim) = false;
                    if any(group)
                        traces_cis = cell2mat(cellfun(@(x) x(neuro_idx,:),Traces(group),'uni',0));
                        mn_sz = min([size(traces,2), size(traces_cis,2)]);
                        dist_cis(irep) = nanmean(reshape(pdist2(traces(:,randperm(rseed,size(traces,2),mn_sz))',...
                            traces_cis(:,randperm(rseed,size(traces_cis,2),mn_sz))'),[],1));
                    else
                        dist_cis(irep) = nan;
                    end
                    
                    dist_mean(irep) = pdist2(sphmean(traces)',sphmean(traces_trans)');
                    smean(irep) = diangle(sphmean(traces));
                    svar(irep) = sphvar(traces);
                    rmean(irep) = mean(sqrt(sum(traces.^2)));
                    rvar(irep) = var(sqrt(sum(traces.^2)));
                end
                
                key.corr = nanmean(corrs);
                key.dist_in = nanmean(dist_in);
                key.dist_trans = nanmean(dist_trans);
                key.dist_cis = nanmean(dist_cis);
                key.dist_mean = nanmean(dist_mean);
                key.smean = nanmean(smean);
                key.svar = nanmean(svar);
                key.rmean = nanmean(rmean);
                key.rvar = nanmean(rvar);
                
                % insert
                key.movie_name = Stims{iStim};
                self.insert(key)
                
            end
        end
    end
    
    methods
        function [Data, Stims, info, Unit_ids] = getData(self,key,bin,stim_split)
            
            if nargin<3 || isempty(bin)
                bin = fetch1(obj.SiteStatsOpt & key, 'binsize');
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
        
        function plot(self,type,varargin)
            
            params.trained = false;
            params.mask = false;
            
            params = getParams(params,varargin);
            
            [data,brain_areas,movie_names, animal_ids,neurons,keys] = fetchn(self,type,'brain_area','movie_name','animal_id','neurons');
            un_areas = unique(brain_areas);
            R = [];
            stim = cellfun(@(x) str2num(x{1}),regexp(movie_names,'\d(?=\w*v\d)','match'));
            is_trained = false(size(stim));
            if params.trained
                for animal = fetch(mice.Mice & self,'animal_id')'
                    tstim = 0;
                    if count(beh.MovieClipCond & animal)
                        tstim = cellfun(@(xx) str2num(cell2mat(xx)),regexp(unique(fetchn(beh.MovieClipCond & animal,'movie_name')),'\d(?=\w*v\d)','match'))';
                    end
                    idx = animal_ids==animal.animal_id;
                    is_trained(idx) = any(stim(idx)'==tstim');
                end
                leg_name = {'Naive','Trained'};
            end
            for iarea = 1:length(un_areas)
                for itrained = 1:max(is_trained)+1
                    R{iarea,itrained} = data(strcmp(un_areas{iarea},brain_areas) & is_trained==itrained-1);
                end
            end
            
            % plot
            figure
            if params.mask
                colors = parula(30);
                plotMask(anatomy.Masks)
                colormap parula
                c = colorbar;

                % set min/max
                mx = max(cellfun(@nanmean,R));
                mn = min(cellfun(@nanmean,R));

                for iarea = 1:length(un_areas)
                    mi = nanmean(R{iarea});
                    if isnan(mi);continue;end
                    idx = double(uint8(floor(((mi-mn)/(mx - mn))*0.99*size(colors,1)))+1);
                    plotMask(anatomy.Masks & ['brain_area="' un_areas{iarea} '"'],colors(idx,:),sum(~isnan(R{iarea})))
                end           
                set(c,'ytick',linspace(0,1,5),'yticklabel',roundall(linspace(mn,mx,5),10^-round(abs(log10(mn/10)))))
                ylabel(c,type,'Rotation',-90,'VerticalAlignment','baseline','fontsize',10)
            else
                b = boxfun(R);
                set(gca,'xtick',1:size(R,1),'xticklabel',un_areas)
                ylabel(type)
                if params.trained && max(is_trained)>0
                    l = legend(b(:,1),leg_name);
                    set(l,'box','off')
                end
            end
        end
        
    end
    
    methods(Static)
        function pwr = fitTune(traces)
            pwr = nan(size(traces,2),1);
            for itrace = 1:size(traces,2)
                x = sort(traces(:,itrace),'descend');
                y = repmat((1:size(x,1))',1,size(x,2));
                [xData, yData] = prepareCurveData( y, x );
                
                % Set up fittype and options.
                ft = fittype( 'power2' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                opts.Display = 'Off';
                opts.StartPoint = [0 0 0];
                
                % Fit model to data.
                [fitresult, ~] = fit( xData, yData, ft, opts );
                v = coeffvalues(fitresult);
                pwr(itrace) = abs(v(2));
            end
            
        end
    end
end
