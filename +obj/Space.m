%{
# Distances in high dimentional space
-> obj.Dec
---
hyp_dist                 : longblob               # distance from the hyperplane
raw_dist                 : longblob               # avg distance between stimuli
stim_idx                 : mediumblob             # stimuli index
%}

classdef Space < dj.Computed
    
    properties
        keySource  = obj.Dec
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % get data
            [classifier,trial_info] = fetch1(obj.Dec & key,'classifier','trial_info');
            [Traces, Stims, ~, Unit_ids] = getData(obj.Dec,key,[],0);
            
            % get distances to hyperplanes
            stim_names = cellfun(@unique,trial_info.names,'uni',0);
            [~,stim_idx] = intersect(Stims,[stim_names{:}]);
            Traces = Traces(stim_idx);
            Stims = Stims(stim_idx);

            % build stim index vector
            nclasses = length(stim_names);
            stimVec = nan(1,length(Stims));
            for iStim = 1:nclasses
                 [~,stim_sub] = intersect(Stims,stim_names{iStim});
                 stimVec(stim_sub) = iStim;
            end
            stimVec = cell2mat(cellfun(@(x,y) ones(1,size(x,2))*y,Traces,num2cell(stimVec),'uni',0));
            Traces = cell2mat(Traces);
                        
            hDist = [];eDist = [];
            for icellnum = 1:size(classifier{1}{1},1)
                for irep = 1:size(classifier,2)
                    
                    % find correct cell indexes
                    [~,~,cell_idxB] = intersect(trial_info.units{irep}{icellnum},Unit_ids,'stable');
                    dat = Traces(cell_idxB,:);
                    
                    % distances in raw space
                    Dis = squareform(pdist(dat'));
                    
                    for iclass = 1:nclasses
                        
                        % distances from hyperplane
                        if iclass==2 && nclasses==2; planeSign=-1;else planeSign=1;end
                        svm = classifier{irep}{iclass}{icellnum};
                        HypDist = @(X) dot((X - svm(3,:)')./svm(4,:)' ,repmat(svm(1,:)',1,size(X,2)),1) + svm(2,1);
                        hDist{iclass}(irep,:,icellnum) = planeSign * single(HypDist(dat));
                        
                        % get mean distances between different classes
                        eDist{iclass}(irep,:,icellnum) = single(nanmean(Dis(stimVec==iclass,:)));
                        
                    end
                end
            end
            
            % insert
            key.hyp_dist = hDist;
            key.raw_dist = eDist;
            key.stim_idx = stimVec;
            self.insert(key)
        end
    end
    
    methods
        
        function [Stims, info] = getInfo(self,key,norepeats)
                        
            if nargin<3|| isempty(norepeats)
                norepeats = true;
            end
            % fetch stimuli without repeats
            if norepeats
                trial_obj = stimulus.Trial &  ...
                    ((stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"')) - ...
                    (aggr(stimulus.Clip , stimulus.Trial & key, 'count(*)->n') & 'n>1')) & key;
            else
                trial_obj = stimulus.Trial &  ...
                    ((stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"'))) & key;
            end
            [flip_times, trial_idxs] = fetchn(...
                trial_obj,'flip_times','trial_idx','ORDER BY trial_idx');
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            trial_idxs = trial_idxs(tidx);
            Stims = unique(fetchn(stimulus.Clip & trial_obj,'movie_name'));
            
            % split for unique stimuli
            for istim = 1:length(Stims)
                [s_trials,s_clips,s_names] = fetchn(stimulus.Trial * stimulus.Clip & ...
                    sprintf('movie_name = "%s"',Stims{istim}) & key, 'trial_idx','clip_number','movie_name');
                [tr_idx, b]= ismember(trial_idxs,s_trials);
                st_idx = b(b>0);
                ti = fetch1(obj.Dec & key,'trial_info');
                info.bins{istim} = reshape(repmat(1:max(ti.bins{1}),sum(tr_idx),1)',[],1);
                info.trials{istim} = reshape(repmat(s_trials(st_idx),1,max(ti.bins{1}))',[],1);
                info.clips{istim} = reshape(repmat(s_clips(st_idx),1,max(ti.bins{1}))',[],1);
                info.names{istim} = reshape(repmat(s_names(st_idx),1,max(ti.bins{1}))',[],1);
            end
        end
        
        function cell_num = getCellNum(self, key)
            if nargin<2; restr = proj(self);else restr = key;end
            ti = fetchn(obj.Dec & restr,'trial_info');
            cell_num = [];
            for ikey = 1:length(ti)
                cell_num{ikey} = cellfun(@length,ti{ikey}.units{1});
            end
            if ikey==1;cell_num=cell_num{1};end
        end
        
        function stim_names = getStim(self, key)
            if nargin<2; restr = proj(self);else restr = key;end
            ti = fetchn(obj.Dec & restr,'trial_info');
            stim_names = [];
            for ikey = 1:length(ti)
                stim_names{ikey} = cellfun(@unique,ti{ikey}.names,'uni',0);
            end
            if ikey==1;stim_names=stim_names{1};end
        end
    end
end
