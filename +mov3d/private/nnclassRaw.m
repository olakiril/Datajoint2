function mi = nnclassRaw(traces,varargin)

% function [CA,CR,FP,FN] = nnclassRaw(traces)
%
% performs a nearest neighbor classification
% and outputs mutual information
% plus the false positives,false negatives,correct acceptance and
% correct rejections.
% traces: [cells classes trials]
%
% MF 2012-11-29

p.avg = 1;

% update parameters if supplied
p = getParams(p,varargin);

% get the sizes
ntrials = size(traces,3);

% initialize
pairs = nchoosek(1:size(traces,2),2);
mi = cell(size(pairs,1),1);
nclasses = 2;

% loop through the pairs
for ipair = 1:size(pairs,1)
    data = traces(:,pairs(ipair,:),:);
    
    % initialize
    F = zeros(nclasses);
    [CA,CR,FP,FN] = initialize('zeros',nclasses,1);
    
    % loop through trials
    for iTrial = 1:ntrials
        
        % calculate mean without taking that trial into account
        ind = true(ntrials,1);
        ind(iTrial) = false;
        r = nanmean(data(:,:,ind),3);
        
        % loop through classes
        for iClass = 1:nclasses
            [foo,indx] = min(pdist2(r',data(:,iClass,iTrial)'));
            F(iClass,indx) = F(iClass,indx) + 1;
        end
    end
    
    % loop through classes
    d = diag(F,0);
    for iclass = 1:nclasses
        CA(iclass) = F(iclass,iclass);
        dind = true(size(d));dind(iclass) = false;
        CR(iclass) = sum(d(dind));
        FN(iclass) = sum(F(iclass,dind));
        FP(iclass) = sum(F(dind,iclass));
    end
    CM = zeros(2,2);
    CM(1,1) = sum(CA);
    CM(1,2) = sum(FN);
    CM(2,1) = sum(FP);
    CM(2,2) = sum(CR);
    
    pA = CM/sum(CM(:));
    pi = sum(CM,2)/sum(CM(:));
    pj = sum(CM,1)/sum(CM(:));
    pij = pi*pj;
    if FN+FP == 0 % this is wrong, it should be FN+FP
        mi{ipair} = 1;
    elseif CA+CR == 0
        mi{ipair} = 0;
    else
        mi{ipair} = sum(sum(pA.*log2(pA./pij)));
    end
end

if p.avg
    mi = nanmean(cell2mat(mi));
else
    mi = cell2mat(mi);
end

%% old method
% % initialize
% F = zeros(nclasses);
% [CA,CR,FP,FN] = initialize('zeros',nclasses,1);
%
% % loop through trials
% for iTrial = 1:ntrials
%
%     % calculate mean without taking that trial into account
%     ind = true(ntrials,1);
%     ind(iTrial) = false;
%     r = mean(traces(:,:,ind),3);
%
%     % loop through classes
%     for iClass = 1:nclasses
%         [~,indx] = min(pdist2(r',traces(:,iClass,iTrial)'));
%         F(iClass,indx) = F(iClass,indx) + 1;
%     end
% end
%
% % loop through classes
% d = diag(F,0);
% for iclass = 1:nclasses
%     CA(iclass) = F(iclass,iclass);
%     dind = true(size(d));dind(iclass) = false;
%     CR(iclass) = sum(d(dind));
%     FN(iclass) = sum(F(iclass,dind));
%     FP(iclass) = sum(F(dind,iclass));
% end
%
% CM(1,1) = sum(CA);
% CM(1,2) = sum(FN);
% CM(2,1) = sum(FP);
% CM(2,2) = sum(CR);
%
% p = CM/sum(CM(:));
% pi = sum(CM,2)/sum(CM(:));
% pj = sum(CM,1)/sum(CM(:));
% pij = pi*pj;
%
% mi = sum(sum(p.*log2(p./pij)));