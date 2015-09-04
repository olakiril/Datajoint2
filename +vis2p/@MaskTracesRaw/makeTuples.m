function makeTuples( obj, key,traces,traces2,coordinates,masknums)

import vis2p.*

% detect & remove duplicates
thr = 5; % threshold in microns
mtraces = mean(traces); % get mean of traces
d = squareform(pdist(coordinates)); % get pairwise distances
d(logical(tril(ones(size(d,1),size(d,2)),0)))=thr; % select upper triangle
[indx,indy] = ind2sub(size(d),find(d<thr));ind=[indx indy];% get indexes
[~,sind]=sort(mtraces(ind),2,'ascend'); % select low firing ones
traces(:,ind(sind==1))      = []; % remove extra cells
traces2(:,ind(sind==1))     = []; % remove extra cells

for itrace = 1:size(traces,2)
    tuple=key;
    tuple.masknum = masknums(itrace);
    try
        mtype = fetch1(MaskTraces(tuple),'mask_type');
    catch
        mtype = 'neuron';
    end
    tuple.mask_type = mtype;
    tuple.calcium_trace = single(traces(:,itrace));
    try
         tuple.red_trace = single(traces2(:,itrace));
    end
    insert( obj, tuple );
end

