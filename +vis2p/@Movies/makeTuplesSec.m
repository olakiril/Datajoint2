function makeTuples( obj, key )

tuple = key;

if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
   %% MPScan 
    % Load movie information
    tpr = tpReader(Scans(key));
    
    % Align movies
    if ~isAligned(tpr)
        %         return
        display('Aligning movie');
        close(tpr);
        fileNameCore = getLocalPath([fetch1(Experiments(key),'directory') '/' fetch1(Scans(key),'file_name')]);
        alignFile_calculate(fileNameCore);
        alignFile_apply(fileNameCore);
        createFilePreview(fileNameCore);
        tpr = tpReader(Scans(key));  % re-open
    end
    
    % do the rest
    channel = getImagingChannel(tpr,1);
    [tuple.xsize, tuple.ysize, tuple.nframes] = size(channel);
    tuple.fps = getFramerate( tpr );
    
    %  Compute mean raw frames
    disp('computing unaligned mean frames');
    tuple.raw_green = calcmean(channel,tuple.nframes);
    channel = getImagingChannel(tpr,2);
    tuple.raw_red = calcmean(channel,tuple.nframes);
    if isAligned( tpr )
        % Compute mean aligned frames
        disp('computing aligned mean frames');
        channel = getAlignedImagingChannel(tpr,1);
        tuple.affine_green = calcmean(channel,tuple.nframes);
        channel = getAlignedImagingChannel(tpr,2);
        tuple.affine_red = calcmean(channel,tuple.nframes);
        
        % get the aligning info
        fileName = [fetch1(Scans(key),'file_name') '_alignment.mat'];
        path = fetch1(Experiments(key),'directory');
        alignData = load(getLocalPath([path '/' fileName]));
        tuple.alignXshifts = alignData.alignData.shifts(1,:);
        tuple.alignYshifts = alignData.alignData.shifts(2,:);
    end
    
    
    % insert patched cell with masknum=0
    if  strcmp(fetch1(Scans(key),'aim'),'patching') && ~isempty(getElectrodeChannelIndicies(tpr) == 2)
        
        % get ephys data
        elChan  = getElectrodeChannel(tpr,2);
        tuple.ephys_fs = getSamplingRate(elChan);
    end
    
    insert( obj, tuple );
    
elseif strcmp(fetch1(Scans(key),'scan_prog'),'ScanImage')
    %% ScanImage file
    
    % Load movie information
    flnm = getFilename(Scans(key));
    scim = tpMethods.Reader(flnm{1});
    
    disp 'reading tiff file'
    movie = scim.read(1);
    gmean = mean(movie,3);
    gmean = gmean-min(gmean(:));
    gmean = uint8(255*gmean./max(gmean(:)));

    disp 'raster correction'
warp = ne7.micro.RasterCorrection.fit(double(g), [3 5]);
tuple.raster_correction = warp;
g = ne7.micro.RasterCorrection.apply(g, warp);

disp 'motion correction...'
assert(scim.hdr.acq.fastScanningX==1 & scim.hdr.acq.fastScanningY==0, 'x must be the fast axis')

offsets = ne7.micro.MotionCorrection.fit(g);
offsets = bsxfun(@minus, offsets, median(offsets));
tuple.motion_correction = int16(offsets);
tuple.motion_rms = mean(std(offsets).*[tuple.um_height tuple.um_width]./[tuple.px_height tuple.px_width]);

disp 'averaging frames...'
g = ne7.micro.MotionCorrection.apply(g, offsets);
tuple.green_img = single(mean(g,3));
    
    % calculate data
    [tuple.xsize, tuple.ysize, tuple.nframes] = size(movie);
    tuple.fps = scim.fps;
    tuple.raw_green = calcmean(movie,tuple.nframes);
    tuple.raw_red = calcmean(movie,tuple.nframes);
    tuple.affine_green = calcmean(cmovie,tuple.nframes);
    tuple.affine_red = calcmean(cmovie,tuple.nframes);
    
    % insert patched cell with masknum=0
    if  strcmp(fetch1(Scans(key),'aim'),'patching') && scim.hasChannel(4)
        tuple.ephys_fs = scim.lps;
    end
    
    insert( obj, tuple );
    
elseif strcmp(fetch1(Scans(key),'scan_prog'),'AOD')
    %% AOD Scan
    [path,name] = fetch1( Experiments*Scans(key), 'directory','file_name' );
    if isempty(Scans(key,'exp_date > "2012-05-01"'))
        filename = getLocalPath([path '/' name '.h5']);
    else
        
        dirs = dir(['M:\Mouse\' name(1:10) '*']);
        names = vertcat(dirs.name);
        % find the session that started immediatelly before the recording
        timediff = str2num(names(:,[12 13 15 16 18 19]))- str2num( name([12 13 15 16 18 19]));
        timeind = find(timediff<0);
        [foo,i] = max(timediff(timeind));
        itime = timeind(i);
        filename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' name];
    end
    % find volumes
    u_ind = strfind(name,'_');
    volumename = fetchn(Scans(rmfield(key,'scan_idx'),[...
        'file_name LIKE "' name(1:u_ind(1)) ...
        '%" and scan_idx > ' num2str(key.scan_idx - 2) ...
        ' and scan_idx < ' num2str(key.scan_idx + 2) ...
        ' and aim = "Stack"' ...
        ]),'file_name');
    if isempty(volumename)
        volumename = fetchn(Scans(rmfield(key,'scan_idx'),[...
            'file_name LIKE "' name(1:u_ind(1)) ...
            '%" and scan_idx > ' num2str(key.scan_idx - 3) ...
            ' and scan_idx < ' num2str(key.scan_idx + 3) ...
            ' and aim = "Stack"' ...
            ]),'file_name');
    end
    if isempty(volumename)
        display('no volumes coupled to the file!,skipping...')
        s = 0;
        x =0;
        y = 0;
        z = 0;
        s2 = 0;
    else
        if ~isempty(Scans(key,'exp_date > "2012-05-01"'))
            volumename = ['M:\Mouse\' dirs(i).name '\AODAcq\' volumename{1}];
            fn = aodReader(volumename,'Volume');
            
            s = fn(:,:,:,1);
            x = fn.x;
            y =  fn.y;
            z =  fn.z;
            s2 = [];
        else
            volumename = getLocalPath([path '/' volumename{1} '.h5']);
            [s x y z s2] = loadAODStack(volumename);
        end
    end
    
    if ~isempty(Scans(key,'exp_date > "2012-05-01"'))
        [traces t] = loadAODTraces(filename);
    else
        [traces t] = loadTraces(filename);
    end
    
    try
        [xpos ypos zpos times] = trackMotion(filename);
        tuple.alignXshifts = interp1(times,xpos,t2);
        tuple.alignYshifts = interp1(times,ypos,t2);
        tuple.alignZshifts = interp1(times,zpos,t2);
        
    catch
        
        tuple.alignXshifts = 0;
        tuple.alignYshifts = 0;
        tuple.alignZshifts = 0;
        
    end
    tuple.xsize        = length(x);
    tuple.ysize        = length(y);
    tuple.zsize        = length(z);
    tuple.nframes      = size(traces,1);
    tuple.fps          = 1/mean(diff(t));
    tuple.raw_green    = s;
    tuple.raw_red      = s2;
    
    insert( obj, tuple );
   
elseif strcmp(fetch1(Scans(key),'scan_prog'),'Unirec')
     %% Unirec Scan
    [path,name] = fetch1( Experiments*Scans(key), 'directory','file_name' );
    filename = getLocalPath([path '/' name '%d.h5']);
    br = baseReader(filename);
    
    tuple.alignXshifts = 0;
    tuple.alignYshifts = 0;
    tuple.alignZshifts = 0;
    tuple.xsize        = 1;
    tuple.ysize        = 1;
    tuple.zsize        = 1;
    tuple.nframes      = getNbSamples(br);
    tuple.fps          = getSamplingRate(br);
    tuple.raw_green    = [];
    tuple.raw_red      = [];
    
    insert( obj, tuple );
end

% calculate mean without excess of memory
function meanCh = calcmean(channel,nframes)
frameBlock = 1000;
meanCh = 0;

% no need for mean of the whole file  %%%% CHECK!
if nframes>1000
    nframes = 1000;
end

for iFrame=1:frameBlock:nframes
    meanCh = meanCh + sum(double(channel(:,:,iFrame:min(iFrame+frameBlock-1,nframes))),3);
end
meanCh = single(meanCh/nframes);