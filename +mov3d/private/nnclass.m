function [J, rJ] = nnclass(traces,varargin)

% function [J rJ] = nnclass(traces)
%
% performs a nearest neighbor classification
% traces: [cells classes trials]
%
% MF 2011-08-25

params.repetitions = 1;
params.cells = size(traces,1);
params.frames = 0;
params.trials = 0;

params = getParams(params,varargin);

if params.frames; y=size(traces,2);else y=1;end

if params.cells == size(traces,1)
    params.repetitions = 1;
end

% do it
J = nan(length(params.cells),y);
rJ = nan(length(params.cells),1);
ic = 0;
for iCell = params.cells;
    p = nan(params.repetitions,size(traces,3),size(traces,2));
    for iRep = 1:params.repetitions
        cellindx = randperm(size(traces,1));
        data = traces(cellindx(1:iCell),:,:);
        for iTrial = 1:size(traces,3)
            ind = true(size(traces,3),1);
            ind(iTrial) = false;
            r = mean(data(:,:,ind),3);
            for iClass = 1:size(traces,2)
                dist = pdist2(r',data(:,iClass,iTrial)');
                [~,indx] = min(dist);
                p(iRep,iTrial,iClass) = indx == iClass;
            end
        end
    end
    ic = ic+1;
    if params.trials
        J = squeeze(p)';
    else
        if params.frames; J(ic,:) = squeeze(mean(mean(p,2),1))';
        else  J(ic) = mean(p(:));
        end
    end
end

if nargout>1
    % randomize traces
    randtraces = traces(:);
    rindx = 1:numel(traces);
    for i = 1:1000
        rindx = rindx(randperm(numel(traces)));
    end
    randtraces = reshape(randtraces(rindx),size(traces,1),size(traces,2),size(traces,3));
    
    rJ = nan(length(params.cells),1);
    ic = 0;
    for iCell = params.cells;
        rp = nan(params.repetitions,size(traces,3),size(traces,2));
        for iRep = 1:params.repetitions
            cellindx = randperm(size(traces,1));
            randdata = randtraces(cellindx(1:iCell),:,:);
            for iTrial = 1:size(traces,3)
                ind = true(size(traces,3),1);
                ind(iTrial) = false;
                rr = mean(randdata(:,:,ind),3);
                for iClass = 1:size(traces,2)
                    dist = pdist2(rr',randdata(:,iClass,iTrial)');
                    [~,rindx] = min(dist);
                    rp(iRep,iTrial,iClass) = rindx == iClass;
                end
            end
        end
        ic = ic+1;
        rJ(ic) = mean(rp(:));
    end
end