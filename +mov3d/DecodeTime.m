%{
mov3d.DecodeTime (computed) # calcium trace
-> preprocess.Sync
-> mov3d.DecodeOpt
-> preprocess.SpikeMethod
-> preprocess.Method
---
obj_cis                      : longblob                      # classification performance for trials that follow same object
obj_trans                    : longblob                      # classification performance for trials that follow opposite object
%}

classdef DecodeTime < dj.Relvar & dj.AutoPopulate
    %#ok<*AGROW>
    %#ok<*INUSL>
    
    properties
        popRel  = (experiment.Scan  ...
            * (preprocess.Spikes & 'spike_method = 3'  & 'extract_method=2'))...
            * (mov3d.DecodeTimeOpt & 'process = "yes"') ...
            * (preprocess.Sync & (vis.MovieClipCond & (vis.Movie & 'movie_class="object3d"')))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            [dec_method,trial_bins,trial_method] = ...
                fetch1(mov3d.DecodeTimeOpt & key,...
                'decode_method','trial_bins','trial_method');
            
            [Data, ~, ~, ~, Trials] = getData(self,key); % [Cells, Obj, Trials]
                       
            % create trial index
            switch trial_method
                case 'random'
                    trial_idx = randperm(size(Data,3));
                case 'sequential'
                    trial_idx = 1:size(Data,3);
            end
            trial_bin = floor(size(Data,3)/trial_bins);
            minBinNum = mode(diff(find(abs(diff(reshape(Trials',[],1)))>0)));
            obj_cis = cell(trial_bins,1);
            obj_trans = cell(trial_bins,1);
            % run the decoding
    
            parfor itrial = 1:trial_bins
                
                display(['Decoding trial # ' num2str(itrial)])
                tIdx = trial_idx(...
                    1+trial_bin*(itrial-1):trial_bin*itrial);
                data = Data(:,:,tIdx);
                trials = Trials(:,tIdx);
              
                if strcmp(dec_method,'nnclassSV')
                   % mi = eval([dec_method '(data,''trials'',1)']);
                   mi = nnclassSV(data,'trials',1);
                end
                utrials = unique(trials(:));
                
                cis = [];
                trans = [];
                for iseg = 2:length(utrials)
                    idx = find(trials==utrials(iseg));
                    [x,~] = ind2sub(size(trials),idx);
                    [xbef,~] = ind2sub(size(trials),find(trials==utrials(iseg-1)));
                    if all(xbef(end)==x)
                        cis{end+1} = mi(idx);
                    elseif all(xbef(end)~=x)
                        trans{end+1} = mi(idx);
                    end
                end
                idx = cellfun(@(x) length(x)==minBinNum,cis);
                obj_cis{itrial} = mean(cell2mat(cis(idx)),2);
                idx = cellfun(@(x) length(x)==minBinNum,trans);
                obj_trans{itrial} = mean(cell2mat(trans(idx)),2);
            end
            
            % insert
            tuple.obj_cis = obj_cis;
            tuple.obj_trans = obj_trans;
            
            % correct for key mismach
            tuple = rmfield(tuple,'spike_method');
            tuple = rmfield(tuple,'extract_method');
            tuple.spike_inference = key.spike_method;
            tuple.segment_method = key.extract_method;
            self.insert(tuple)
        end
    end
    
    methods
        function [Data, xloc, yloc, zloc, Trials] = getData(obj,key,ibin)
         
            bin = fetch1(mov3d.DecodeTimeOpt & key, 'binsize');
            if nargin>2;bin = ibin;end
            
          
            
            AA = []; BB = [];
            
            nslices = fetch1(preprocess.PrepareGalvo & key, 'nslices');
            
            for islice = 1:nslices
                key.slice = islice;
                
                caTimes = fetch1(preprocess.Sync &  (experiment.Scan & key), 'frame_times');
                
                caTimes = caTimes(key.slice:nslices:end);
                [X, xloc{islice}, yloc{islice}, zloc{islice}] = ...
                    fetchn(preprocess.SpikesRateTrace * preprocess.MaskCoordinates & key,...
                    'rate_trace','xloc','yloc','zloc');
                X = [X{:}];
                xm = min([length(caTimes) length(X)]);
                X = @(t) interp1(caTimes(1:xm)-caTimes(1), X(1:xm,:), t, 'linear', nan);  % traces indexed by time
                
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
                AA{islice} = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),A,'uni',0)'),size(A{1},2),[]),[3 1 2]);
                AA_trials{islice} = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),A_trials,'uni',0)'),size(A_trials{1},2),[]),[3 1 2]);
                
                B = snippet(2,:);
                B_trials = trial_idx(2,:);
                idx = ~cellfun(@isempty,B);
                B = B(idx);
                B_trials = B_trials(idx);
                BB{islice} = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),B,'uni',0)'),size(B{1},2),[]),[3 1 2]);
                BB_trials{islice} = permute(reshape(cell2mat(cellfun(@(x) reshape(x',[],1),B_trials,'uni',0)'),size(B_trials{1},2),[]),[3 1 2]);
            end
            
            % Arrange data
            objA = cell2mat(AA);
            objA_trials = cell2mat(AA_trials);
            objB = cell2mat(BB);
            objB_trials = cell2mat(BB_trials);
            mS = min([size(objA,3) size(objB,3)]);
            Data = reshape(permute(objA(:,:,1:mS),[2 4 3 1]),size(objA,2),1,[]);
            Data(:,2,:) = reshape(permute(objB(:,:,1:mS),[2 4 3 1]),size(objB,2),1,[]);
            Trials = reshape(permute(objA_trials(:,:,1:mS),[2 4 3 1]),size(objA_trials,2),1,[]);
            Trials(:,2,:) = reshape(permute(objB_trials(:,:,1:mS),[2 4 3 1]),size(objB_trials,2),1,[]);
            Trials = squeeze(Trials(1,:,:));
        end
   
        function plot(obj,k)
            bin = fetch1(mov3d.DecodeTimeOpt & k,'binsize');
            obj_cis = fetch1(mov3d.DecodeTime & k,'obj_cis');
            obj_trans = fetch1(mov3d.DecodeTime & k,'obj_trans');
            idx = ~cellfun(@isempty,obj_cis);
            cis = cell2mat(obj_cis(idx)');
            idx = ~cellfun(@isempty,obj_trans);
            trans = cell2mat(obj_trans(idx)');
            figure
            errorPlot(1:size(cis,1),cis','errorColor',[0 0 0.5'])
            hold on
            errorPlot(1:size(cis,1),trans','errorColor',[0.5 0 0])
            ylim([0.5 1])
            set(gca, 'xtick',1:size(cis,1),'xticklabel',bin/1000:bin/1000:5) 
        end
    end
    
end
