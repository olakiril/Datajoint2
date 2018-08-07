%{
# Movie repeats evaluation
-> obj.RepeatsOpt
-> fuse.ScanDone
---
-> stimulus.Sync
%}

classdef Repeats < dj.Computed
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        keySource = fuse.ScanDone * (obj.RepeatsOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial & (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"'))))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            % get parameters & data
            disp 'fetching data ...'
            [method, shuffle] = fetch1(obj.RepeatsOpt & key,'method','shuffle');
            [Data, Hashes, keys] = getData(self,key);
            
            % insert
            tuple = key;
            self.insert(tuple)
            
            % initialize
            disp 'computing reliabity ...'
            import_keys = cell(size(Data,4),1);
            rl_sh = cell(size(Data,4),1);
            rl = cell(size(Data,4),1);
            parfor icell = 1:size(Data,4)
                % shuffling
                cell_traces = Data(:,:,:,icell); %[trials time stim]
                rand_idx = cell2mat(arrayfun(@(x) randperm(numel(cell_traces)),1:shuffle,'uni',0)');
                r_cell_traces = reshape(cell_traces(rand_idx),[shuffle size(cell_traces)]);
                clip_keys = [];
                
                for iclip = 1:size(cell_traces,3)
                    trace = cell_traces(:,:,iclip); %[trials time]
                    switch method
                        case 'explainedVar'
                            rl{icell} = var(nanmean(trace,2))/nanvar(trace(:));
                            rl_sh{icell} = var(nanmean(r_cell_traces(:,:,:,iclip),3),[],2)/nanvar(trace(:));
                        case 'corr'
                            traceZ = zscore(trace);
                            r_traceZ = permute(zscore(r_cell_traces(:,:,:,iclip),[],2),[2 3 1]);
                            c = corr(traceZ);
                            sel_idx = logical(tril(ones(size(traceZ,2)),-1));
                            rl{icell} = nanmean(c(sel_idx));
                            rh_sh{icell} = nan(shuffle,1);
                            for i = 1:shuffle
                                c = corr(r_traceZ(:,:,i));
                                rh_sh{icell}(i) = nanmean(c(sel_idx));
                            end
                        case 'oracle'
                            traceZ = zscore(trace);
                            r_traceZ = zscore(r_cell_traces(:,:,:,iclip),[],2);
                            c = nan(size(traceZ,2),1);r_c = nan(size(traceZ,2),shuffle);
                            for irep = 1:size(traceZ,2)
                                idx = true(size(traceZ,2),1);
                                idx(irep) = false;
                                c(irep) = corr(traceZ(:,irep),mean(traceZ(:,idx),2));
                                r_c(irep,:) = diag(corr(r_traceZ(:,:,irep)',mean(r_traceZ(:,:,idx),3)'));
                            end
                            rl{icell} = nanmean(c);
                            rl_sh{icell} = nanmean(r_c,2);
                    end
                    
                    % clear tuple
                    tuple = keys(icell);
                    tuple.rep_opt = key.rep_opt;
                    tuple.condition_hash = Hashes(iclip).condition_hash;
                    tuple.r = rl{icell};
                    tuple.r_shuffle = mean(rl_sh{icell});
                    tuple.p_shuffle = mean(rl_sh{icell}>rl{icell});
                    clip_keys{iclip} = tuple;
                end
                import_keys{icell} = cell2mat(clip_keys);
            end
            
            % insert
            disp 'Inserting ...'
            for tuple = cell2mat(import_keys')
                insert(obj.RepeatsUnit,tuple)
            end
        end
    end
    
    methods
        function [Data,Hashes,keys] = getData(self,key,bin) % [time,repeats,uni_stims,cells]
            
            if nargin<3 || isempty(bin)
                bin = fetch1(obj.RepeatsOpt & key, 'binsize');
            end
            
            % get stuff
            [Traces, caTimes, keys] = getAdjustedSpikes(fuse.ActivityTrace & key,'soma'); % [time cells]
            Traces = single(Traces);
            trials = stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d"')) & (obj.Repeats & key);
            [flip_times, cond_hash] = fetchn(...
                trials * (stimulus.Clip & (aggr(stimulus.Clip, stimulus.Trial & (obj.Repeats & key), 'count(*)->n') & 'n>1')),...
                'flip_times','condition_hash');
            
            % filter out incomplete trials
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            
            % find repeated trials
            [uni_hashes,~,uni_idx] = unique(cond_hash);
            
            % Subsample traces
            flip_times = cell2mat(flip_times(tidx));
            xm = min([length(caTimes) length(Traces)]);
            X = @(t) interp1(caTimes(1:xm)-caTimes(1), Traces(1:xm,:), t, 'linear', nan);  % traces indexed by time
            fps = 1/median(diff(flip_times(1,:)));
            d = max(1,round(bin/1000*fps));
            traces = convn(permute(X(flip_times - caTimes(1)),[2 3 1]),ones(d,1)/d,'same');
            traces = traces(1:d:end,:,:);
            
            % split for unique stimuli
            Data = [];Hashes = struct('condition_hash',[]);
            for istim = unique(uni_idx)'
                Hashes(istim).condition_hash = uni_hashes{istim};
                Data{istim} = permute(traces(:,:,uni_idx==istim),[1 4 3 2]);
            end
            
            % equilize trials
            Data = permute(cell2mat(cellfun(@(x) x(:,:,1:min(cellfun(@(x) size(x,3),Data)),:),Data,'uni',0)),[1 3 2 4]);
            
        end
        
        function plot(self,key,bin,varargin)
            
            params.contrast = 0.5;
            params.gap = 50;
            params.scroll_step = 3;
            params.trials = 0;
            params.rand = 0;
            params.sort = 0;
            params.prctile = 0;
            
            params = getParams(params,varargin);
            
            if nargin<3 || isempty(bin)
                bin = fetch1(self.RepeatsOpt & key, 'binsize');
            end
            
            if nargin<2 || isempty(key)
                key = fetch(self);
%             else
%                 keys = fetch(self & key);
             end
%             for key = keys'
                    [Data, Stims] = getData(self,key,bin); % [time,repeats,uni_stims,cells]
                    dat= []; s = []; im = [];
                    fig = figure;
                    for iobj = 1:size(Data,3)
                        data = permute(squeeze(Data(:,:,iobj,:)),[3,1,2]); % [Cells, Time, Trials]
                        if params.sort
                            [~,mxidx] = max(nanmean(data,3),[],2);
                            [~,sidx] = sort(mxidx,'ascend');
                            data = data(sidx,:,:);
                        end
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
%             end
            
            function plotData
                clf
                set(fig,'WindowScrollWheelFcn', @scroll,'KeyPressFcn',@dispkeyevent)
                s = [];
                for iobj = 1:length(dat)
                    s(iobj) = subplot(1,length(dat),iobj);
                    im(iobj) = imagesc(abs(dat{iobj}'.^params.contrast));
                    
                    colormap (1 - gray)
                    step_per_cell = trial_num(iobj)+floor(params.gap/100*trial_num(iobj));
                    set(gca,'ytick',trial_num(iobj)/2:step_per_cell:dat_sz(iobj),...
                        'yticklabel',1:size(data,1),...
                        'box','off',...
                        'xtick',0:1000/bin: bin_num(iobj),...
                        'xticklabel',0:1:bin/1000* bin_num(iobj))
                    [mov,clip] = fetch1(stimulus.Clip & Stims(iobj),'movie_name','clip_number');
                    title(sprintf('Movie: %s\nclip: %d',mov,clip));
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
