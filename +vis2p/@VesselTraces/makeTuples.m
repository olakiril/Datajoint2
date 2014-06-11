function makeTuples( obj, key)

% get scan information
[lens,mag] = fetch1(Scans(key),'lens','mag');
pixelPitch = 11000/512/mag/lens;

% get data
tpr = tpReader(Scans(key));
if ~ isAligned(tpr);display('Not Aligned! Skipping...');return;end
gChan = getAlignedImagingChannel( tpr, 1 );
times = readSyncedTimestamps(tpr);
data = gChan(:,:,:);

% invert movies & downsample by 10X
data = -data + max(data(:));
idx = 1:10:size(data,3)-10;
for iFrame =  1:length(idx)
    data(:,:,iFrame) = mean(data(:,:,idx(iFrame) : idx(iFrame)+10),3);
    times(iFrame) = times(idx(iFrame));
end
data = data(:,:,1:length(idx));
times = times(1:length(idx));

% find the good vessels across all frames
im = mean(data,3);
cellFinder = CellFinder( im, pixelPitch, 'minRadius', 2, 'minDistance', 20,...
    'minSharpness',1,'minContrast', 0.05,'maxCells',30,'nodisplay',1);
[x1,y1,contrast,sharpness,correlation] = getCells( cellFinder );

% find the vessels for each frame
[x,y,radii] = initialize('nan',size(data,3),length(x1));
for iFrame = 1:size(data,3)
    if ~mod(iFrame,100)
        display([num2str(iFrame) '/' num2str(size(data,3))])
    end
    im = data(:,:,iFrame);
    cellFinder = CellFinder( im, pixelPitch, 'minRadius', 2, 'minDistance', 20,...
        'minSharpness',1,'minCorr',0.7,'minContrast', 0.001,'maxCells',40,'nodisplay',1);
    [xt,yt,radiit] = getCells( cellFinder );
    [xi, yi] = ind2sub([length(xt) length(x1)],find(pdist2([xt yt],[x1 y1])<10));
    x(iFrame,yi) = xt(xi);
    y(iFrame,yi) = yt(xi);
    radii(iFrame,yi) = radiit(xi);
end

% insert data
for iCell = 1:size(x,2);
    tuple=key;
    tuple.vmasknum = iCell;
    tuple.x_trace = single(x(:,iCell));
    tuple.y_trace = single(y(:,iCell));
    tuple.radii_trace = single(radii(:,iCell));
    tuple.contrast= contrast(iCell);
    tuple.sharpness = sharpness(iCell);
    tuple.fit_qual = correlation(iCell);  % correlation between the fitted blob and the pixel values
    tuple.timestamps = times;
    insert( obj, tuple );
end

