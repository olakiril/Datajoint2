function makeTuples( obj, key )

[classes classrep bin type repetitions] = fetch1(StatsModelParams(key),'classes','classrep','bin','type','repetitions');
traces = getTraces(StatsSites(key,'binsize = 100'),'bin',bin);

if isempty(traces)
    display('Nothing to compute!')
    return
end

% measure mean and variance of the responses
m = mean(traces,3);
v = var(traces,[],3);

switch type
    case 'rand' % randomize traces
        rindx = 1:numel(traces);
        for i = 1:1000
            rindx = rindx(randperm(numel(traces)));
        end
        traces = reshape(traces(rindx),size(traces,1),size(traces,2),size(traces,3));
        
    case 'cov' % get rid of covariance structure by randomizing trials
        traces = randdim(traces,3);
        
    case 'var' % randomize variance across mean
        rindx = 1:numel(traces);
        for i = 1:1000
            rindx = rindx(randperm(numel(traces)));
        end
        pm = repmat(m,[1 1 size(traces,3)]);
        vr = traces(:) - pm(:);
        traces = pm + reshape(vr(rindx),size(traces,1),size(traces,2),size(traces,3));
end

% select classes
if size(traces,2)<classes
    classes = size(traces,2);
end

% do it
J = nan(size(traces,1),1);
cP = nan(size(traces,1),classes);
parfor iCell = 1:size(traces,1);
    p = nan(repetitions,classrep,size(traces,3),classes);
    for iClassRep = 1:classrep
        Data = traces(:,randperm(size(traces,2),classes),:); % randomly select classes (when classes<max(classes))
        for iRep = 1:repetitions
            data = Data(randperm(size(traces,1),iCell),:,:); %randomly select a group of cells
            for iTrial = 1:size(traces,3)
                r = mean(data(:,:,1:size(traces,3)~=iTrial),3); %mean population response, excluding the test trial
                for iClass = 1:classes
                    dist = pdist2(r',data(:,iClass,iTrial)');
                    [~, indx] = min(dist);
                    p(iRep,iClassRep,iTrial,iClass) = indx == iClass;
                end
            end
        end
    end
    J(iCell) = mean(p(:));
    cP(iCell,:) = mean(mean(mean(permute(p,[4 3 2 1]),4),3),2);
end
    
key.class_mean = m;
key.class_var = v;
key.performance = J;
key.class_per = cP(end,:);

% insert data
insert( obj, key );
