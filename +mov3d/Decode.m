%{
mov3d.Decode (computed) # calcium trace
-> rf.Sync
-> mov3d.DecodeOpt
-> pre.SpikeInference
-> pre.SegmentMethod
---
mi                    : longblob                      # mutual information
%}

classdef Decode < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = (rf.Scan & pre.Spikes) * (pre.SpikeInference & 'spike_inference = 2')...
            * pre.SegmentMethod * (mov3d.DecodeOpt & 'process = "yes"') * ...
            (rf.Sync & psy.Movie & (psy.Session & 'stimulus="object3d"') )
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            AA = [];
            BB = [];
            tuple = key;
            
            nslices = fetch1(pre.ScanInfo & key, 'nslices');
            [bin,sel_method,dec_method,trial_bins,trial_method] = fetch1(mov3d.DecodeOpt & key,...
                'binsize','select_method','decode_method','trial_bins','trial_method');
            for islice = 1:nslices
                key.slice = islice;
                
                caTimes = fetch1(rf.Sync &  (rf.Scan & key), 'frame_times');
                
                caTimes = caTimes(key.slice:nslices:end);
                [X,xloc{islice},yloc{islice},zloc{islice}] = fetchn(pre.Spikes * pre.MaskCoordinates & key,...
                    'spike_trace','xloc','yloc','zloc');
                X = [X{:}];
                xm = min([length(caTimes) length(X)]);
                X = @(t) interp1(caTimes(1:xm)-caTimes(1), X(1:xm,:), t, 'linear', nan);  % traces indexed by time
                
                trials = pro(rf.Sync*psy.Trial & (rf.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
                trials = fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx');
                
                snippet = [];
                stims = [2 1];
                for trial = trials'
                    stim = stims(~isempty(strfind(trial.path_template,'obj1'))+1);
                    % extract relevant trace & bin
                    fps = 1/median(diff(trial.flip_times));
                    t = trial.flip_times - caTimes(1);
                    d = max(1,round(bin/1000*fps));
                    trace = convn(X(t),ones(d,1)/d,'same');
                    snippet{stim,end+1} = trace(1:d:end,:);
                end
                
                A = snippet(1,:);
                A = A(~cellfun(@isempty,A));
                AA{islice} = reshape(cell2mat(A),size(A{1},1),size(A{1},2),[]);
                
                B = snippet(2,:);
                B = B(~cellfun(@isempty,B));
                BB{islice} = reshape(cell2mat(B),size(B{1},1),size(B{1},2),[]);
            end
            
            % Arrange data
            objA = cell2mat(AA);
            objB = cell2mat(BB);
            mS = min([size(objA,3) size(objB,3)]);
            Data = reshape(permute(objA(:,:,1:mS),[2 4 3 1]),size(objA,2),1,[]);
            Data(:,2,:) = reshape(permute(objB(:,:,1:mS),[2 4 3 1]),size(objB,2),1,[]);
            
            % compute distances for 'expand' method
            xloc = cell2mat(xloc');yloc = cell2mat(yloc');zloc = cell2mat(zloc');
            
            % create trial index
            switch trial_method
                case 'random'
                    trial_idx = randperm(size(Data,3));
                case 'sequential'
                    trial_idx = 1:size(Data,3);
            end
            trial_bin = floor(size(Data,3)/trial_bins);
            if strcmp(sel_method,'all');mi = nan(trial_bins,1);else mi = nan(trial_bins,size(Data,1));end
            
            % run the decoding
            for itrial = 1:trial_bins
                data = Data(:,:,trial_idx(...
                    1+trial_bin*(itrial-1):trial_bin*itrial));
                
                switch sel_method
                    case 'all'
                        mi(itrial) = eval([dec_method '(data)']);
                    case 'subsample'
                        cell_idx = randperm(size(Data,1));
                        for icell = 1:size(Data,1)
                            dat = data(cell_idx(1:icell),:,:); %#ok<NASGU>
                            mi(itrial,icell) =  eval([dec_method '(dat)']);
                        end
                    case 'expand'
                        cell_idx = randperm(size(Data,1),1);
                        [~,sort_idx] = sort(abs(pdist2([xloc(cell_idx),yloc(cell_idx),zloc(cell_idx)],[xloc,yloc,zloc])),'ascend');
                        dat = data(sort_idx,:,:); %#ok<NASGU>
                        for icell = 1:size(Data,1)
                            mi(itrial,icell) =  eval([dec_method '(dat(1:icell,:,:))']);
                        end
                end
            end
            
            % insert
            tuple.mi = mi;
            self.insert(tuple)
            
        end
    end
    
    methods
        function Data = getData(obj,key,bin)
            
            if nargin<3;bin = 500;end
            
            AA = [];
            BB = [];
            
            nslices = fetch1(pre.ScanInfo & key, 'nslices');
            
            for islice = 1:nslices
                key.slice = islice;
                
                caTimes = fetch1(rf.Sync &  (rf.Scan & key), 'frame_times');
                
                caTimes = caTimes(key.slice:nslices:end);
                [X,xloc{islice},yloc{islice},zloc{islice}] = fetchn(pre.Spikes * pre.MaskCoordinates & key,...
                    'spike_trace','xloc','yloc','zloc');
                X = [X{:}];
                xm = min([length(caTimes) length(X)]);
                X = @(t) interp1(caTimes(1:xm)-caTimes(1), X(1:xm,:), t, 'linear', nan);  % traces indexed by time
                
                trials = pro(rf.Sync*psy.Trial & (rf.Scan & key) & 'trial_idx between first_trial and last_trial', 'cond_idx', 'flip_times');
                trials = fetch(trials*psy.Movie, '*', 'ORDER BY trial_idx');
                
                snippet = [];
                stims = [2 1];
                for trial = trials'
                    stim = stims(~isempty(strfind(trial.path_template,'obj1'))+1);
                    % extract relevant trace & bin
                    fps = 1/median(diff(trial.flip_times));
                    t = trial.flip_times - caTimes(1);
                    d = max(1,round(bin/1000*fps));
                    trace = convn(X(t),ones(d,1)/d,'same');
                    snippet{stim,end+1} = trace(1:d:end,:);
                end
                
                A = snippet(1,:);
                A = A(~cellfun(@isempty,A));
                AA{islice} = reshape(cell2mat(A),size(A{1},1),size(A{1},2),[]);
                
                B = snippet(2,:);
                B = B(~cellfun(@isempty,B));
                BB{islice} = reshape(cell2mat(B),size(B{1},1),size(B{1},2),[]);
            end
            
            % Arrange data
            objA = cell2mat(AA);
            objB = cell2mat(BB);
            mS = min([size(objA,3) size(objB,3)]);
            Data = reshape(permute(objA(:,:,1:mS),[2 4 3 1]),size(objA,2),1,[]);
            Data(:,2,:) = reshape(permute(objB(:,:,1:mS),[2 4 3 1]),size(objB,2),1,[]);
        end
        
        function plotMasks(obj,key)
            figure
            colors = parula(50);
            plotMask(map.Masks)
            colormap parula
            colorbar
            
            
            areas =  fetchn(map.Area,'area');
            for iarea = 1:length(areas)
                
                keys = fetch(obj & key & (rf.Scan & ['cortical_area="' areas{iarea} '"']));
                if isempty(keys);continue;end
                for ikey = 1:length(keys)
                    tuple = keys(ikey);
                    
                    mi(ikey) = mean(fetch1(obj & tuple,'mi'));
                    
                end
                mi = mean(mi);
                %                  idx = ceil(mi*50);
                idx = ceil(mi*100)-50;
                plotMask(map.Masks & ['area="' areas{iarea} '"'],colors(idx,:))
            end
            
        end
        
        function plotMasksNorm(obj)
            key.dec_opt = 1;
            sessions = fetch(rf.Session & obj);
            areas =  fetchn(map.Area,'area');
            
            MI = nan(length(areas),length(sessions));
            
            for isession = 1:length(sessions)
                for iarea = 1:length(areas)
                    keys = fetch(obj & key & sessions(isession) & (rf.Scan & ['cortical_area="' areas{iarea} '"']));
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
            
            figure
            colors = parula(100);
            plotMask(map.Masks)
            colormap parula
            c = colorbar;
            for iarea = 1:length(areas)
                if isnan(nmi(iarea));continue;end
                idx = ceil(((nmi(iarea) - mn_nmi) /(mx_nmi  - mn_nmi))*100);
                if idx ==0;idx =1;end
                plotMask(map.Masks & ['area="' areas{iarea} '"'],colors(idx,:))
            end
            set(gcf,'name','Normalized performance')
            c.Ticks = 0:0.2:1;
            c.TickLabels = roundall(mn_nmi: (mx_nmi - mn_nmi)/5 :mx_nmi,0.1);
            ylabel(c, 'Normalized performance to V1')
            
        end
        
        function plot(obj,key,colors,linestyle)
            if nargin<2; key = [];end
            if ~isfield(key,'dec_opt')
                key.dec_opt = 11;
            end
            
            if ~isfield(key,'spike_inference')
                key.spike_inference = 2;
            end
            
            if ~isfield(key,'segment_method')
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
            for idx = 1:length(keys)
                tuple = keys(idx);
                mi = fetch1(mov3d.Decode & tuple,'mi');
                [name,name2] = fetch1(rf.Scan & tuple,'cortical_area','scan_notes');
                if strcmp(name,'other');name = name2;end
                errorPlot(1:size(mi,2),mi,'errorColor',colors(idx,:),'linestyle',linestyle);
                names{idx} = name;
            end
            
            
            xlabel('Neuron #')
            ylabel('Mutual Information (bits)')
            set(gca,'box','off')
%             l = legend(names);
%             set(l,'box','off','location','northwest')
            title('Classifier: SVM, Bin size = 0.5sec')
        end
        
        function plotLLE(obj,k)
%             k.animal_id = 9508;
%             k.spike_inference = 2;
%             k.segment_method = 2;
%             k.session = 1;
%             k.scan_idx = 2;
            
            traces = getData(mov3d.Decode,k,5000);
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