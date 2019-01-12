%{
# Distances in high dimentional space
-> obj.Dec
---
score                    : mediumblob             # distance from the hyperplane
distance                 : mediumblob             # avg distance between stimuli
stim_idx                 : mediumblob             # stimuli index
stim_class               : blob                   # stimuli classes
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
            key.score = hDist;
            key.distance = eDist;
            key.stim_idx = stimVec;
            key.stim_class = Stims;
            self.insert(key)
        end
    end
    
end
