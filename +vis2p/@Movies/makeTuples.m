function makeTuples( obj, key )

import vis2p.*

tuple = key;

if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
    %% MPScan
    % Load movie information
    tpr = tpReader(Scans(key));
    fileNameCore = getLocalPath([fetch1(Experiments(key),'directory') '/' fetch1(Scans(key),'file_name')]);

    % Align movies
    if ~exist(getLocalPath([fileNameCore '_fineAlignment.mat']),'file') && false
        display('Aligning movie');
        im = getImagingChannel(tpr,1);
        mc.warp_degree = 4;
        mc.xwarp_degree = 0;
        yWarp = tpMethods.YWarp( double(mean(im(:,:,1:100),3)));
        degrees = [mc.warp_degree*[1 1] mc.xwarp_degree];
        mc.warp_polynom = zeros(size(im,3), sum(degrees)+2, 'single');
        p = zeros(1, sum(degrees)+2);
        tic
        for iFrame = 1:size(im,3)
            if ~mod(iFrame,100);fprintf('[%3d/%d]\n',iFrame,size(im,3));end
            % fit polynomials
            yWarp.fit(double(im(:,:,iFrame)), degrees, p);
            p = yWarp.coefs;
            mc.warp_polynom(iFrame, :) = p;
        end
        save(getLocalPath([fileNameCore '_fineAlignment.mat']),'mc')
        toc

    elseif ~isAligned(tpr) && ~exist([fileNameCore '_alignment.mat'],'file')
        display('Aligning movie');
        close(tpr);
        fileNameCore = getLocalPath([fetch1(Experiments(key),'directory') '/' fetch1(Scans(key),'file_name')]);
        alignFile_calculate(fileNameCore);
        try
            alignFile_apply(fileNameCore);
            createFilePreview(fileNameCore);
        catch
            display('Cannot save file!')
        end
        tpr = tpReader(Scans(key));  % re-open
    end
    
    % do the rest
    channel1 = getImagingChannel(tpr,1);
    [tuple.xsize, tuple.ysize, tuple.nframes] = size(channel1);
    tuple.fps = getFramerate( tpr );
    
    %  Compute mean raw frames
    disp('computing unaligned mean frames');
    tuple.raw_green = calcmean(channel1,tuple.nframes);
    channel2 = getImagingChannel(tpr,2);
    tuple.raw_red = calcmean(channel2,tuple.nframes);
    
    disp('computing aligned mean frames');
    if exist(getLocalPath([fileNameCore '_fineAlignment.mat']),'file')
        % get the aligning info
        load(getLocalPath([fileNameCore '_fineAlignment.mat']));
        tuple.motion_correction = mc;
                
        % apply alignment
        degrees = [mc.warp_degree mc.warp_degree 0];
        frames = zeros(size(channel1,1),size(channel1,2),length(1:2:size(channel1,3)),'uint16');
        for i=1:2:size(channel1,3)
            frames(:,:,i) = uint16(tpMethods.YWarp.apply(double(channel1(:,:,i)), mc.warp_polynom(i,:), degrees));
        end
        tuple.affine_green = calcmean(frames,tuple.nframes); 

        frames = zeros(size(channel2,1),size(channel2,2),length(1:2:size(channel2,3)),'uint16');
        for i=1:2:size(channel2,3)
            frames(:,:,i) = uint16(tpMethods.YWarp.apply(double(channel2(:,:,i)), mc.warp_polynom(i,:), degrees));
        end
        tuple.affine_red = calcmean(frames,tuple.nframes); 
        
    elseif isAligned( tpr )
        % Compute mean aligned frames 
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
    elseif exist([fileNameCore '_alignment.mat'],'file') && ~isAligned(tpr)
        % get the aligning info
        fileName = [fetch1(Scans(key),'file_name') '_alignment.mat'];
        path = fetch1(Experiments(key),'directory');
        alignData = load(getLocalPath([path '/' fileName]));
        tuple.alignXshifts = alignData.alignData.shifts(1,:);
        tuple.alignYshifts = alignData.alignData.shifts(2,:);
        
        % apply alignment
        stack = getData(channel1);
        stack = transformStack(stack, alignData.alignData.transMatrices);
        validRect = alignData.alignData.validRect;
        stack = stack(validRect(3):validRect(4), validRect(1):validRect(2), :);
        tuple.affine_green = calcmean(stack,tuple.nframes);
    
        stack = getData(channel2);
        stack = transformStack(stack, alignData.alignData.transMatrices);
        validRect = alignData.alignData.validRect;
        stack = stack(validRect(3):validRect(4), validRect(1):validRect(2), :);
        tuple.affine_red = calcmean(stack,tuple.nframes);
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
    
    % get mouse info to see whether to compute channel 2
    transMice = {'NesCre-ReYFP','Viat-Ai9','SST-Ai9','PV-Ai9','PV-AAVArch','SST-AAVArch','Nestin-Ai9'};
    rDyes = {'FL4,SR','FL4,CR','SR','CR','OGB,SR','OGB,SR,CR','OGB,CR'};
    redMouse = any(strcmp(fetch1(Mice(key),'mouse_strain'),transMice));
    redChannel = any(strcmp(fetch1(Experiments(key),'dyes'),rDyes));
    
    % Load movie information
    scim = tpReader(Scans(key));
    
    disp 'reading tiff file'
    % check for channel 2
    if scim.hasChannel(2) && redMouse || scim.hasChannel(2) && redChannel
        movie = double(scim.read([1 2]));
        redmovie = movie(:,:,:,2);
        movie = movie(:,:,:,1);
    else
        movie = double(scim.read(1));
    end
    [tuple.xsize, tuple.ysize, tuple.nframes] = size(movie);
    tuple.raw_green = calcmean(movie,tuple.nframes);
    tuple.fps = scim.fps;
    
    disp 'raster correction'
    raster = tpMethods.RasterCorrection.fit(movie, [3 5]);
    tuple.raster_correction = raster;
    rmovie = single(tpMethods.RasterCorrection.apply(movie, raster));
    clear movie 
    
    disp 'motion correction...'
    assert(scim.hdr.acq.fastScanningX==1 & scim.hdr.acq.fastScanningY==0, 'x must be the fast axis')
    offsets = tpMethods.MotionCorrection.fit(rmovie);
    offsets = bsxfun(@minus, offsets, median(offsets));
    movie = tpMethods.MotionCorrection.apply(rmovie, offsets);
    tuple.motion_correction = offsets;
    
    disp 'averaging frames...'
    tuple.affine_green = single(mean(movie,3));
    
%         disp 'computing subpixel correction...'
%         tuple.motion_correction.warp_degree = 2;
%         tuple.motion_correction.xwarp_degree = 4;
%         yWarp = tpMethods.YWarp(tuple.affine_green);
%         degrees = [tuple.motion_correction.warp_degree*[1 1]...
%             tuple.motion_correction.xwarp_degree];
%         tuple.motion_correction.warp_polynom = ...
%             zeros(size(movie,3), sum(degrees)+2, 'single');
%         p = zeros(1, sum(degrees)+2);
%         for iFrame = 1:size(movie,3)
%             if ~mod(sqrt(iFrame),1);fprintf('[%3d/%d]\n',iFrame,size(movie,3));end
%     
%             % fit polynomials
%             yWarp.fit(double(rmovie(:,:,iFrame)), degrees, p);
%             p = yWarp.coefs;
%             tuple.motion_correction.warp_polynom(iFrame, :) = p;
%         end
%         tuple.affine_green = calcmean(getFrames(Movies, 1, [], ...
%             tuple.raster_correction, tuple.motion_correction,movie),nframes);
    
    % do it for second channel
    if scim.hasChannel(2) && redMouse || scim.hasChannel(2) && redChannel
        disp 'red channel...'
        movie = redmovie;
        tuple.raw_red = calcmean(movie,tuple.nframes);
        movie = tpMethods.RasterCorrection.apply(movie, raster);
        movie = tpMethods.MotionCorrection.apply(movie, offsets);
        tuple.affine_red = single(mean(movie,3));
    else
        tuple.raw_red = [];
        tuple.affine_red = [];
    end
    
    % insert patched cell with masknum=0
    if  strcmp(fetch1(Scans(key),'aim'),'patching') && scim.hasChannel(4)
        tuple.ephys_fs = scim.lps;
    end
    
    % indert
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
            timediff = str2num(names(:,[12 13 15 16 18 19]))- str2num( volumename{1}([12 13 15 16 18 19]));
            timeind = find(timediff<0);
            [~,i] = max(timediff(timeind));
            itime = timeind(i);
        
            volumename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' volumename{1}];
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