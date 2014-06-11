function makeTuples( obj, key )

% get possitions
[cellx, celly, cellz] = fetchn(MaskCells(key),'img_x','img_y','img_z');

% get Traces
[traces,~,T] = getTraces(StatsSites(rmfield(key,'neurons')).*StatsSitesParams(key,'stats_opt = 1'));
masknums = fetchn(T,'masknum');
traces = traces(:,:)';

% binarize
s = std(traces);
traces = bsxfun(@(t,s) t > key.thr_factor * s,traces,s);

%% run analysis
euDist = @(x1,x2,y1,y2,z1,z2) sqrt((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2);

    
nrns = size(traces,2);

% manipulate distances
if strcmp(key.dist_fact,'random')
    seed = key2num(rmfield(key,{'movie_type','thr_factor','edge_crop', ...
        'trace_opt','binsize','dist_fact'}));
    rand('state',seed); %#ok<RAND>
    indx = randperm(size(cellx,1));
    cx = cellx(indx);
    cy = celly(indx);
    cz = cellz(indx);
else
    f = str2double(key.dist_fact);
    cx = cellx.*exp(f);
    cy = celly.*exp(f);
    cz = cellz.*exp(-f);
end

for iNeuron = 1:nrns
    
    tuple = key;
    tuple.masknum = masknums(iNeuron);
    trueDist = euDist(cellx(iNeuron),cellx,celly(iNeuron),celly,cellz(iNeuron),cellz);
    cellDist = euDist(cx(iNeuron),cx,cy(iNeuron),cy,cz(iNeuron),cz);
    [dist, sortseq] = sort(cellDist); %#ok<ASGLU>
    indx = sortseq(1:key.neurons);
    
    tuple.r = runPopulationAnalysis(traces(:,indx)');
    tuple.mdist = mean(trueDist(indx(2:end)));
    
    % insert data
    insert( obj, tuple );
end


