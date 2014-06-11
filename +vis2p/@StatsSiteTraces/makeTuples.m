function makeTuples( obj, key )

% get traces of all neurons of a site and fps
T = Traces(key,'masknum > 0');

% select cells
k.mouse_id = key.mouse_id;
k.exp_date= key.exp_date;
k.scan_idx= key.scan_idx;
if ~isempty(MaskTracesQuality(k))
    T = T.*MaskTracesQuality(k,['ca_event_snr>' num2str(key.trace_qual)]);
end
k.trace_opt= key.trace_opt;
if ~isempty(RFFit(k))
    T = T.*RFFit(k,['snr >' num2str(key.rf_snr_thr)]);
end
if ~isempty(RFStats(k))
    T = T.*RFStats(k,['onpoff_p <' num2str(key.rf_p_thr)]);
end
tracesR = fetchn(T,'trace');

% get
tracesR = [tracesR{:}];
fps = fetch1( Movies(key), 'fps' );

% get positions 
[cellx celly cellz] = fetchn(MaskCells.*T,'img_x','img_y','img_z');

% Load times of the traces
times = fetchn(VisStims(key),'frame_timestamps');
times = reshape(cell2mat(times),1,[]);

% stim times and types of movies shown
stimTimes = fetchn(StatsPresents(key),'movie_times');
movieTypes = fetchn(StatsPresents(key),'movie_num');
movieTypes = bsxfun(@eq,movieTypes,unique(movieTypes)');

% find trace segments for each stimulus
traces = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = tracesR(times > (stimTimes{iTimes}(1)) & ...
        times < (stimTimes{iTimes}(end)),:);
end

if isempty(traces)
    display('Nothing to compute!')
    return
end

% remove incomplete trials
L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
indx = L(:,1) >= mean(L(:,1))*9/10;
traces = traces(indx);
L = L(indx);
movieTypes = movieTypes(indx,:);

% equalize
traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);

% rebin to approximate binsize (downsampling by an integer factor)
traces = cat(3,traces{:});
d = max(1,round(key.binsize/1000*fps));
key.actual_binsize = d*1000/fps;
k = ones(d,1)/d;
traces = convn(traces,k,'valid');
traces = double(traces(1:d:end,:,:));

% mean across same stimulus of different repetitions and collapse segments
t = cell(2,1);
s = cell(2,1);
for iMovie = 1:size(movieTypes,2)
    t{iMovie} = mean(traces(:,:,movieTypes(:,iMovie)),3);
    z = zscore(traces(:,:,movieTypes(:,iMovie)),[],3);
    s{iMovie} = reshape(permute(z,[1 3 2]),[],size(z,2));
end
mtraces = cell2mat(t);
straces = cell2mat(s);
clear t
clear s

%index build
[iCell,jCell] = ndgrid(1:size(mtraces,2),1:size(mtraces,2));

% compute distances
euDist = @(x1,x2,y1,y2,z1,z2) sqrt((x1 - x2).^2 + (y1 - y2).^2 + (z1 - z2).^2);
cellIndx = jCell(iCell>jCell);
cellIndy = iCell(iCell>jCell);
key.site_traces_dist = euDist(cellx(cellIndx),cellx(cellIndy),celly(cellIndx),...
    celly(cellIndy),cellz(cellIndx),cellz(cellIndy));

% compute signal correlogram
c = corr(mtraces);
c = c(iCell>jCell);
key.signal_cor = c ;

% compute noize correlogram
nc = corr(straces);
nc = nc(iCell>jCell);
key.noise_cor = nc ;

% compute angle
key.site_angle = nan(length(c),1);
for i = 1:length(c)
    key.site_angle(i) = acos(dot(mtraces(:,cellIndx(i)),mtraces(:,cellIndy(i))) / ...
        (norm(mtraces(:,cellIndx(i)))*norm(mtraces(:,cellIndy(i)))));
end

% insert data
insert( obj, key );

