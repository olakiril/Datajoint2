%{
# Movie repeats evaluation
-> obj.RepeatsOpt
-> fuse.ScanDone
---
-> stimulus.Sync
r                    : longblob                      # reliability
r_shuffle            : longblob                      # chance
p_shuffle            : longblob                      # p value
%}

classdef Repeats < dj.Computed
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        keySource  = fuse.ScanDone * (obj.RepeatsOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"'))))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            [method, shuffle] = fetch1(obj.RepeatsOpt & key,'method','shuffle');
            
            [Data,~,rData] = getData(self,key,[],shuffle); % [Cells, Time, Trials]
            
            switch method
                case 'explainedVar'
                    r = cell2mat(cellfun(@reliability,Data, mat2cell(ones(size(Data)),1,ones(size(Data))),'uni',0));
                    parfor ishuffle = 1:shuffle
                        r_shuffle(:,:,ishuffle) = cell2mat(cellfun(@reliability,cellfun(@(x) x(:,:,:,ishuffle),rData,'uni',0), mat2cell(ones(size(rData)),1,ones(size(rData))),'uni',0));
                    end
                    p_shuffle = nan(size(r));
                    for i = 1:size(r,1)
                        for y = 1:size(r,2)
                            p_shuffle(i,y) = mean(r(i,y)<r_shuffle(i,y,:));
                        end
                    end
                case 'corr'
                    r = nan(size(Data{1},1),length(Data));
                    for iobj = 1:length(Data)
                        cor = nan(size(Data{iobj},1),1);
                        for icell = 1:size(Data{iobj},1)
                            traceZ = zscore(squeeze(Data{iobj}(icell,:,:)));
                            c = corr(traceZ);
                            cor(icell) = nanmean(c(logical(tril(ones(size(traceZ,2)),-1))));
                        end
                        r(:,iobj) = (cor);
                    end
                    r_shuffle = nan(size(Data{1},1),length(Data),shuffle);
                    parfor ishuffle = 1:shuffle
                        r_s = nan(size(rData{1},1),length(rData));
                        for iobj = 1:length(rData)
                            cor = nan(size(rData{iobj},1),1);
                            for icell = 1:size(rData{iobj},1)
                                traceZ = zscore(squeeze(rData{iobj}(icell,:,:,ishuffle)));
                                c = corr(traceZ);
                                cor(icell) = nanmean(c(logical(tril(ones(size(traceZ,2)),-1))));
                            end
                            r_s(:,iobj) = (cor);
                        end
                        r_shuffle(:,:,ishuffle) = r_s;
                    end
                    p_shuffle = nan(size(r));
                    for i = 1:size(r,1)
                        for y = 1:size(r,2)
                            p_shuffle(i,y) = mean(r(i,y)<r_shuffle(i,y,:));
                        end
                    end
                    
                case 'oracle'
                    r = nan(size(Data{1},1),length(Data));
                    for iobj = 1:length(Data)
                        cor = nan(size(Data{iobj},1),1);
                        for icell = 1:size(Data{iobj},1)
                            traceZ = zscore(squeeze(Data{iobj}(icell,:,:)));
                            for irep = 1:size(traceZ,2)
                                idx = true(size(traceZ,2),1);
                                idx(irep) = false;
                                c(irep) = corr(traceZ(:,irep),mean(traceZ(:,idx),2));
                            end
                            cor(icell) = nanmean(c);
                        end
                        r(:,iobj) = (cor);
                    end
                    r_shuffle = nan(size(Data{1},1),length(Data),shuffle);
                    parfor ishuffle = 1:shuffle
                        r_s = nan(size(rData{1},1),length(rData));
                        for iobj = 1:length(rData)
                            cor = nan(size(rData{iobj},1),1);
                            for icell = 1:size(rData{iobj},1)
                                traceZ = zscore(squeeze(rData{iobj}(icell,:,:,ishuffle)));
                                c = nan(size(traceZ,2),1);
                                for irep = 1:size(traceZ,2)
                                    idx = true(size(traceZ,2),1);
                                    idx(irep) = false;
                                    c(irep) = corr(traceZ(:,irep),mean(traceZ(:,idx),2));
                                end
                            cor(icell) = nanmean(c);
                            end
                            r_s(:,iobj) = (cor);
                        end
                        r_shuffle(:,:,ishuffle) = r_s;
                    end
                    p_shuffle = nan(size(r));
                    for i = 1:size(r,1)
                        for y = 1:size(r,2)
                            p_shuffle(i,y) = mean(r(i,y)<r_shuffle(i,y,:));
                        end
                    end
            end
            % insert
            tuple.r = r;
            tuple.r_shuffle = r_shuffle;
            tuple.p_shuffle = p_shuffle;
            self.insert(tuple)
            
        end
    end
    
    methods
        function [Data,Stims,rData] = getData(obj,key,bin,shuffle) % {uni_stims}(cells,time,repeats)
            
            if nargin<3 || isempty(bin)
                [bin] = fetch1(obj.RepeatsOpt & key, 'binsize');
            end
            
            % get stuff
            [Traces, caTimes] = getAdjustedSpikes(fuse.ActivityTrace & key,'soma'); % [time cells]
            trials = stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"')) & key;
            [flip_times, trial_idxs, mov, clip] = (fetchn(trials * stimulus.Clip,'flip_times','trial_idx','movie_name','clip_number'));
            
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
            
            % randomize
            if nargin>3 && shuffle
                sz_data = size(traces);
                rtraces = repmat(reshape(permute(traces,[2 1 3]),sz_data(2),[]),1,1,shuffle);
                rtraces(:,:,1) = rtraces(:,randperm(size(rtraces,2)),1);
                for i = 2:shuffle
                    rtraces(:,:,i) = rtraces(:,randperm(size(rtraces,2)),i-1);
                end
                rtraces = permute(reshape(rtraces,size(rtraces,1),sz_data(1),sz_data(3),shuffle),[2 1 3 4]);
            else
                rtraces = [];
            end
            
            % split for unique stimuli
            Data = [];k=[];Stims = [];rData = [];
            for istim = 1:length(unis)
                Stims{istim,1} = unimov{unis(istim,1)};
                Stims{istim,2} = unis(istim,2);
                k.movie_name = unimov{unis(istim,1)};
                k.clip_number = unis(istim,2);
                stim_trials = fetchn(trials*stimulus.Cond & k,'trial_idx');
                Data{istim} = permute(traces(:,:,ismember(trial_idxs(tidx & ridx),stim_trials)),[2 1 3]);
                if ~isempty(rtraces)
                    rData{istim} = permute(rtraces(:,:,ismember(trial_idxs(tidx & ridx),stim_trials),:),[2 1 3 4]);
                end
            end
        end
        
        function plot(obj,key,bin,varargin)
            
            params.contrast = 0.5;
            params.gap = 50;
            params.scroll_step = 3;
            params.trials = 0;
            params.rand = 0;
            
            params = getParams(params,varargin);
            
            if nargin<3 || isempty(bin)
                bin = fetch1(obj.RepeatsOpt & key, 'binsize');
            end
            
            if nargin<2 || isempty(key)
               keys = fetch(obj);
            else
                keys = fetch(obj.Repeats & key);
            end
            for key = keys'
                try
                    [Data, Stims] = getData(obj,key,bin); % [Cells, Time, Trials]
                    dat= []; s = []; im = [];
                    fig = figure;
                    for iobj = 1:length(Data)
                        data = Data{iobj};
                        if params.rand; data = data(randperm(size(data,1)),:,:);end
                        if params.trials; data = permute(data,[3 2 1]);end
                        trial_num(iobj) = size(data,3);
                        bin_num(iobj) = size(data,2);
                        data(:,:,end+1:end+floor(params.gap/100*trial_num(iobj))) = min(data(:));
                        dat{iobj} = (reshape(permute(data,[2 3 1]),size(data,2),[]));
                        dat_sz(iobj) = size(dat{iobj},2);
                    end
                    step = floor(params.scroll_step/100*min(dat_sz));
                    plotData
                end
            end
            
            function plotData
                clf
                set(fig,'WindowScrollWheelFcn', @scroll,'KeyPressFcn',@dispkeyevent)
                s = [];
                for iobj = 1:length(Data)
                    s(iobj) = subplot(1,length(Data),iobj);
                    im(iobj) = imagesc(dat{iobj}'.^params.contrast);
                    
                    colormap (1 - gray)
                    step_per_cell = trial_num(iobj)+floor(params.gap/100*trial_num(iobj));
                    set(gca,'ytick',trial_num(iobj)/2:step_per_cell:dat_sz(iobj),...
                        'yticklabel',1:size(data,1),...
                        'box','off',...
                        'xtick',0:1000/bin: bin_num(iobj),...
                        'xticklabel',0:1:bin/1000* bin_num(iobj))
                    
                    title(sprintf('Movie: %s\nclip: %d',Stims{iobj,1},Stims{iobj,2}));
                    if iobj ==1
                        xlabel('Time (sec)')
                        if params.trials
                            ylabel('Cells/Trial #')
                        else
                            ylabel('Trials/Cell #')
                        end
                    end
                end
                linkaxes(s,'xy')
                set(gca,'ylim',[0 step_per_cell*20])
            end
            
            function dispkeyevent(~, event)
                switch event.Key
                    case 'f' % toggle fullscreen
                        set(fig,'units','normalized')
                        p = get(fig,'outerposition');
                        if all(p~=[0 0 1 1])
                            set(fig,'outerposition',[0 0 1 1]);
                        else
                            set(fig,'outerposition',f_pos);
                        end
                        set(fig,'units','pixels')
                    case 'comma' % decrease contrast
                        if params.contrast>0.1
                            params.contrast = params.contrast-0.05;
                            for iobj = 1:length(dat)
                                c = get(s(iobj),'children');
                                c.CData = dat{iobj}'.^params.contrast;
                            end
                        end
                    case 'period' % decrease contrast
                        if params.contrast<2
                            params.contrast = params.contrast+0.05;
                            for iobj = 1:length(dat)
                                c = get(s(iobj),'children');
                                c.CData = dat{iobj}'.^params.contrast;
                            end
                        end
                    case 'downarrow'
                        scroll_down
                    case 'uparrow'
                        scroll_up
                    case 'leftarrow'
                        ylim = get(gca,'ylim')*0.9;
                        ylim(ylim<0) = 0;
                        ylim(ylim>min(dat_sz)) = min(dat_sz);
                        set(gca,'ylim',ylim);
                    case 'rightarrow'
                        ylim = get(gca,'ylim')*1.1;
                        ylim(ylim<0) = 0;
                        ylim(ylim>min(dat_sz)) =  min(dat_sz);
                        set(gca,'ylim',ylim);
                    otherwise
                        disp '---'
                        disp(['key: "' event.Key '" not assigned!'])
                end
            end
            
            function scroll(varargin)
                if varargin{2}.VerticalScrollCount<0
                    scroll_up
                elseif varargin{2}.VerticalScrollCount>0
                    scroll_down
                end
            end
            
            function scroll_up
                ylim = get(gca,'ylim');
                if ylim(1) > 1
                    ylim = ylim - min([step ylim(1)]);
                end
                set(gca,'ylim',ylim)
            end
            
            function scroll_down
                ylim = get(gca,'ylim');
                if ylim(2)<min(dat_sz)
                    ylim = ylim + min([step min(dat_sz)-ylim(2)]);
                end
                set(gca,'ylim',ylim)
            end
            
        end
    end
end
