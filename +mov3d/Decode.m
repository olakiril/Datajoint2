%{
mov3d.Decode (computed) # calcium trace
-> preprocess.Sync
-> mov3d.DecodeOpt
-> preprocess.SpikeMethod
-> preprocess.Method
---
mi                    : longblob                      # mutual information
%}

classdef Decode < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  = (experiment.Scan  ...
            * (preprocess.Spikes & 'spike_method = 5'  & 'extract_method=2'))...
            * (mov3d.DecodeOpt & 'process = "yes"') ...
            * (preprocess.Sync  & (vis.MovieClipCond * vis.Trial & ...
            (vis.Movie & 'movie_class="object3d"') & ...
            'trial_idx between first_trial and last_trial'))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            [sel_method,dec_method,trial_bins,trial_method,chance] = ...
                fetch1(mov3d.DecodeOpt & key,...
                'select_method','decode_method','trial_bins','trial_method','chance');
            
            [Data, xloc, yloc, zloc] = getData(self,key); % [Cells, Obj, Trials]
            
            % create trial index
            switch trial_method
                case 'random'
                    trial_idx = randperm(size(Data,3));
                case 'sequential'
                    trial_idx = 1:size(Data,3);
            end
            
            % randomize object identity for chance
            if chance
                sz = size(Data);
                Data = Data(:,:);
                idx = 1:size(Data,2);
                for i = 1:3000
                    idx = idx(randperm(size(Data,2)));
                end
                Data = reshape(Data(:,idx),sz);
            end
            
            
            trial_bin = floor(size(Data,3)/trial_bins);
            if strcmp(sel_method,'all');mi = nan(trial_bins,1);else mi = nan(trial_bins,size(Data,1));end
            
            % run the decoding
            mi = cell(trial_bins,1);
            parfor itrial = 1:trial_bins
                display(['Decoding trial # ' num2str(itrial)])
                
                data = Data(:,:,trial_idx(...
                    1+trial_bin*(itrial-1):trial_bin*itrial));
                
                switch sel_method
                    case 'all'
                        mi{itrial} = decode(data,dec_method);
                    case 'subsample'
                        cell_idx = randperm(size(Data,1));
                        for icell = 1:size(Data,1)
                            dat = data(cell_idx(1:icell),:,:);
                            mi{itrial}(icell) = decode(dat,dec_method);
                        end
                    case 'expand'
                        cell_idx = randperm(size(Data,1),1);
                        [~,sort_idx] = sort(abs(pdist2([xloc(cell_idx),yloc(cell_idx),zloc(cell_idx)],[xloc,yloc,zloc])),'ascend');
                        dat = data(sort_idx,:,:);
                        for icell = 1:size(Data,1)
                            mi{itrial}(icell) = decode(dat(1:icell,:,:),dec_method);
                        end
                end
            end
            
            % insert
            tuple.mi = cell2mat(mi); % [trials cells]
            
            % correct for key mismach
            self.insert(tuple)
        end
    end
    
    methods
        function [Data, xloc, yloc, zloc, Trials, Params] = getData(obj,key,ibin)
            
            [bin, rf_opt] = fetch1(mov3d.DecodeOpt & key, 'binsize','rf_opt');
            if nargin>2;bin = ibin;end
            
            if rf_opt > 0
                index = true;
                rf_key = key;
                rf_key.rf_opt = rf_opt;
                [rf_idx, rf_trials] = fetch1(mov3d.RFFilter & key,'rf_idx','rf_trials');
            else
                index = false;
            end
            
            [Traces, caTimes] = pipetools.getAdjustedSpikes(key);
            [xloc, yloc, zloc] = ...
                fetchn(preprocess.MaskCoordinates & key,...
                'xloc','yloc','zloc');
            xm = min([length(caTimes) length(Traces)]);
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
                if any(isnan(trace(:))) || any(~isreal(trace(:)));continue;end
                if index; trace = trace(rf_idx{trial.trial_idx==rf_trials},:);end
                snippet{stim,idx} = trace;
                trial_idx{stim,idx} = repmat(trial.trial_idx,size(trace,1),1);
            end
            
            A = snippet(1,:);
            Aidx = ~cellfun(@isempty,A);
            A = A(Aidx);
            objA = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),A,'uni',0)'),size(A{1},2),[]),[3 1 2]);
            B = snippet(2,:);
            Bidx = ~cellfun(@isempty,B);
            B = B(Bidx);
            objB = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),B,'uni',0)'),size(B{1},2),[]),[3 1 2]);
            
            % Arrange data
            mS = min([size(objA,3) size(objB,3)]);
            Data = reshape(permute(objA(:,:,1:mS),[2 4 3 1]),size(objA,2),1,[]);
            Data(:,2,:) = reshape(permute(objB(:,:,1:mS),[2 4 3 1]),size(objB,2),1,[]);
            
            % get Trials
            A_trials = trial_idx(1,:);
            A_trials = A_trials(Aidx);
            objA_trials = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),A_trials,'uni',0)'),size(A_trials{1},2),[]),[3 1 2]);
            B_trials = trial_idx(2,:);
            B_trials = B_trials(Bidx);
            objB_trials = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),B_trials,'uni',0)'),size(B_trials{1},2),[]),[3 1 2]);
            Trials = reshape(permute(objA_trials(:,:,1:mS),[2 4 3 1]),size(objA_trials,2),1,[]);
            Trials(:,2,:) = reshape(permute(objB_trials(:,:,1:mS),[2 4 3 1]),size(objB_trials,2),1,[]);
            Trials = squeeze(Trials(1,:,:));
            
            if nargout>5
                % get params
                [params, param_trials] = getParams(obj,key,bin);
                objA_params = cell2mat(params(ismember(param_trials,cellfun(@(x) x(1),A_trials)))');
                objB_params = cell2mat(params(ismember(param_trials,cellfun(@(x) x(1),B_trials)))');
                Params = permute(objA_params(1:mS,:),[2 3 1]);
                Params(:,2,:) = permute(objB_params(1:mS,:),[2 3 1]);
            end
            
        end
        
        function [Data, Stims, info] = getDataNew(obj,key,bin)
            
            if nargin<3
                bin = fetch1(mov3d.DecodeMultiOpt & key, 'binsize');
            end
            
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
            
            % split for unique stimuli
            for istim = 1:length(Stims)
                [s_trials,s_clips,s_names] = fetchn(trials*vis.MovieClipCond &...
                    sprintf('movie_name = "%s"',Stims{istim}),'trial_idx','clip_number','movie_name');
                [tr_idx, b]= ismember(trial_idxs,s_trials);
                st_idx = b(b>0);
                dat = permute(traces(:,:,tr_idx),[2 3 1]);
                info.bins{istim} = reshape(repmat(1:size(dat,3),size(dat,2),1),1,[]);
                info.trials{istim} = reshape(repmat(s_trials(st_idx),1,size(dat,3)),1,[]);
                info.clips{istim} = reshape(repmat(s_clips(st_idx),1,size(dat,3)),1,[]);
                info.names{istim} = reshape(repmat(s_names(st_idx),1,size(dat,3)),1,[]);
                Data{istim} = reshape(dat,size(traces,2),[]);
            end
        end
        
        function [params, param_trials] = getParams(obj,key,bin)
            
            speed = @(x,y,timestep) sqrt(x.^2+y.^2)./timestep;
            binsize= fetch1(mov3d.DecodeOpt & key, 'binsize');
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
        
        function f = plotMasks(obj,norm)
            
            % get data
            method = fetch1(mov3d.DecodeOpt & obj,'decode_method');
            areas =  fetchn(map.Area,'area');
            areas = areas(~strcmp(areas,'MAP'));
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                keys = fetch(obj & (experiment.Scan & ['brain_area="' areas{iarea} '"']));
                if isempty(keys);continue;end
                for ikey = 1:length(keys)
                    tuple = keys(ikey);
                    MI{iarea}(ikey) = mean(fetch1(obj & tuple,'mi'));
                end
            end
            
            % plot
            f = figure;
            colors = parula(30);
            plotMask(map.Masks)
            colormap parula
            c = colorbar;
            if strcmp(method,'nnclassRawSV') || strcmp(method,'nnclassRaw')
                name = 'Mutual information (bits)';
            else
                name = 'Classification performance (%)';
            end
            ylabel(c,name,'Rotation',-90,'VerticalAlignment','baseline')
            
            mxMI = max(cellfun(@nanmean,MI));
            if nargin>1
                mx = mxMI;
            else
                mx =1;
            end
            for iarea = 1:length(areas)
                try
                    mi = nanmean(MI{iarea});
                    if strcmp(method,'nnclassRawSV')
                        idx = ceil(size(colors,1)*0.99/mxMI*mi);
                    else
                        idx = double(uint8(floor(((mi-0.5)/(mx - 0.5))*0.99*size(colors,1)))+1);
                    end
                    plotMask(map.Masks & ['area="' areas{iarea} '"'],colors(idx,:),length(MI{iarea}))
                end
            end
            if strcmp(method,'nnclassSV') || strcmp(method,'nnclass')
                set(c,'ytick',linspace(0,1,5),'yticklabel',roundall(linspace(0.5,mx,5),0.01))
            else
                set(c,'ytick',linspace(0,1,11),'yticklabel',roundall(linspace(0,mxMI,11),0.01))
            end
            
        end
        
        function f = plotMasksNorm(obj)
            sessions = fetch(experiment.Session & obj);
            areas =  fetchn(map.Area,'area');
            areas = areas(~strcmp(areas,'MAP'));
            MI = nan(length(areas),length(sessions));
            
            for isession = 1:length(sessions)
                for iarea = 1:length(areas)
                    
                    keys = fetch(obj & sessions(isession) & (experiment.Scan & ['brain_area="' areas{iarea} '"']));
                    if isempty(keys);continue;end
                    mi = [];
                    for ikey = 1:length(keys)
                        tuple = keys(ikey);
                        mi(ikey) = mean(fetch1(obj & tuple,'mi'));
                    end
                    MI(iarea,isession) = mean(mi);
                end
            end
            nmi = nanmean(bsxfun(@rdivide,MI,MI(1,:)),2);
            mx_nmi = max(nmi);
            mn_nmi = min(nmi);
            
            f = figure;
            colors = parula(100);
            plotMask(map.Masks)
            colormap parula
            c = colorbar;
            for iarea = 1:length(areas)
                try
                    if isnan(nmi(iarea));continue;end
                    idx = ceil(((nmi(iarea) - mn_nmi) /(mx_nmi  - mn_nmi))*100);
                    if idx ==0;idx =1;end
                    plotMask(map.Masks & ['area="' areas{iarea} '"'],colors(idx,:),sum(~isnan(MI(iarea,:))))
                end
            end
            set(gcf,'name','Normalized performance')
            labels = mn_nmi: (mx_nmi - mn_nmi)/5 :mx_nmi;
            c.Ticks = linspace(0,1,length(labels));
            c.TickLabels = roundall(labels,0.01);
            ylabel(c, 'Normalized performance to V1','Rotation',-90,'VerticalAlignment','baseline')
            
        end
        
        function plot(obj,key,colors,linestyle)
            if nargin<2; key = fetch(obj);end
            if ~isfield(key,'dec_opt')
                key.dec_opt = 11;
            end
            
            if ~isfield(key,'spike_method')
                key.spike_inference = 3;
            end
            
            if ~isfield(key,'extract_method')
                key.segment_method = 2;
            end
            
            keys = fetch(mov3d.Decode & key);
            %             figure;
            hold on
            if nargin<3
                colors = hsv(length(keys));
            end
            
            if nargin<4
                linestyle = '-';
            end
            
            names = [];
            figure
            hold on
            for idx = 1:length(keys)
                tuple = keys(idx);
                mi = fetch1(mov3d.Decode & tuple,'mi');
                [name,name2] = fetch1(experiment.Scan & tuple,'brain_area','scan_notes');
                if strcmp(name,'other');name = name2;end
                if size(mi,2)>1
                    errorPlot(1:size(mi,2),mi,'errorColor',colors(idx,:),'linestyle',linestyle);
                else
                    errorbar(idx,mean(mi),std(mi)/sqrt(length(mi)),'.');
                end
                names{idx} = name;
            end
            
            if size(mi,2)>1
                xlabel('Neuron #')
                l = legend(names);
                set(l,'box','off','location','northwest')
            else
                ylim([0 max(get(gca,'ylim'))])
                set(gca,'xtick',1:length(keys),'xticklabel',names)
                xlim([0 length(keys)+1])
            end
            
            ylabel('Mutual Information (bits)')
            set(gca,'box','off')
            title('Classifier: SVM, Bin size = 0.5sec')
        end
        
        function plotLLE(obj,k)
            % bin = 5000;
            bin = 1000;
            if nargin<2
                k.animal_id = 9508;
                k.spike_inference = 2;
                k.segment_method = 2;
                k.session = 1;
                k.scan_idx = 2;
            end
            
            traces = getData(mov3d.Decode,k,bin);
            traces(isnan(traces))=0;
            tracesn = permute(traces,[1 3 2]);
            trac = tracesn(:,:);
            for i = 1:size(trac,1)
                trac(i,:) = conv(trac(i,:),gausswin(10),'same');
            end
            
            Y = lle(trac,12,4);
            
            Y1 = Y(:,1:size(tracesn,2));
            Y2 = Y(:,size(tracesn,2)+1:end);
            
            % Y1 = Y1(:,1:300);
            % Y2 = Y2(:,1:300);
            
            tracs = [];
            for i = 1:size(Y1,1)
                tracs(i,:) = interp1(Y1(i,:),1:0.1:size(Y1,2),'V5cubic');
            end
            Y1 = tracs;
            
            tracs = [];
            
            for i = 1:size(Y1,1)
                tracs(i,:) = interp1(Y2(i,:),1:0.1:size(Y2,2),'V5cubic');
            end
            Y2 = tracs;
            Y1 = Y1([3 1 2],:);
            Y2 = Y2([3 1 2],:);
            
            % plot
            figure
            % PlotTool.sizedFigure([80 40])
            lim = 3;
            ip = size(Y1,2);
            p1 = plot3(Y1(1,1:ip),Y1(2,1:ip),Y1(3,1:ip),'b','linewidth',1);
            grid on
            hold on
            p2 = plot3(Y2(1,1:ip),Y2(2,1:ip),Y2(3,1:ip),'r','linewidth',1);
            
            
            title('LLE V1')
            set(gcf,'name','LLE V1')
        end
    end
    
end
