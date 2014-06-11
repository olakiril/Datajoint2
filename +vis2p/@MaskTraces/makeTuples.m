function makeTuples( obj, key, cellFinder, type,shareMask)

%% ScanImage
if strcmp(fetch1(Scans(key),'scan_prog'),'ScanImage')
    % get site data
    dyes = fetch1(Experiments(key),'dyes');
    [greenImg,redImg] = fetch1( Movies(key),'affine_green','affine_red');
    if isempty(redImg)
        redChannel = false;
    else
        redChannel = true;
        imChan2 = getFrames(Movies(key),2);
    end
    % extract traces
    imChan = getFrames(Movies(key),1);
    
    if nargin>4 && shareMask
        [x,y,z] = fetch1(Scans(key),'x','y','z');
        [xAll,yAll,zAll] = fetchn(Scans(rmfield(key,'scan_idx')).*MaskTraces,'x','y','z');
        idx = find(sum(bsxfun(@eq,[xAll yAll zAll],[x y z]),2)>2);
        if ~isempty(idx)
            keys = fetch(Scans(rmfield(key,'scan_idx')).*MaskTraces);
            shareMask = keys(idx);
        else
            shareMask = false;
        end
    end
    
    if nargin>2 && ~isempty(cellFinder)% manualy insert data
        [cellTraces,~,annulusTraces] = getTraces( cellFinder, imChan );
        redTraces = getTraces( cellFinder, imChan2 );
        
        tuple=key;
        tuple.calcium_trace = single(cellTraces);
        tuple.annulus_trace = single(annulusTraces);
        tuple.red_trace     = single(redTraces);
        if nargin>3
            tuple.mask_type = type;
        else
            tuple.mask_type = 'neuron';
        end
        insert( obj, tuple );
        
        tuple=key;
        [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder ); %#ok<ASGLU>
        greenContrast = getCellContrast( cellFinder, greenImg );
        redContrast   = getCellContrast( cellFinder,   redImg );
        tuple.green_contrast= greenContrast;
        tuple.red_contrast  = redContrast;
        tuple.img_x = x;
        tuple.img_y = y;
        tuple.cell_radius = radii;
        tuple.sharpness = sharpness;
        tuple.fit_qual = correlation;  % correlation between the fitted blob and the pixel values
        makeTuples(MaskCells,tuple)   
    else
        % fit the green image as a sum of radially symmetric blobs
        % with radii proportional to the total magnification
        % find cells
        [lens,mag] = fetch1(Scans(key),'lens','mag');
        disp('Fix this')
        pitch = 11000/mag/lens;
        pixelPitch = pitch/size(imChan,2);
        
        if nargin>4 && isstruct(shareMask) % use masks from another site
            [x,y,radi,contr,sharp,cor] = fetchn(MaskCells(shareMask),...
                'img_x','img_y','cell_radius','green_contrast','sharpness','fit_qual');
            
            cellFinder = CellFinder( greenImg, pixelPitch,...
                'x', x, 'y', y, 'radii', radi,'contrast',contr,...
                'sharpness',sharp,'corrs',cor);
            display(['Using ' shareMask.exp_date ' ' num2str(shareMask.scan_idx) ' Mask'])
        else
            
            cellFinder = CellFinder( greenImg, pixelPitch, 'minRadius', 3.0, 'minDistance', 7.0, 'minContrast', 0.05 );
            fig = figure(1000);
            plot(cellFinder)
            print(fig,'-dpng', [fetch1(Experiments(key),'directory') '\'...
                fetch1(Scans(key),'file_name') '_fcells'])
            
        end
        
        [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder ); %#ok<ASGLU>
        greenContrast = getCellContrast( cellFinder, greenImg );
        [cellTraces,neuropilTrace,annulusTraces] = getTraces( cellFinder, imChan );
        
        if redChannel
            redContrast   = getCellContrast( cellFinder,   redImg );
            isNeuron  =  2.5 > (redContrast-trimmean(redContrast,50))/iqr(redContrast)/0.7413;  % rule of thumb for neurons
            [redTraces,neuropilrTrace] = getTraces( cellFinder, imChan2 );
        else
            isNeuron = true(length(x),1);
        end
        
        % insert neuropil entry with masknum=-1
        tuple=key;
        tuple.masknum=-1;
        tuple.mask_type = 'neuropil';
        tuple.calcium_trace = single(neuropilTrace);
        if redChannel;tuple.red_trace = neuropilrTrace;end
        insert( obj, tuple );
        
        % insert cells entry with masknum=-2
        tuple=key;
        tuple.masknum=-2;
        tuple.mask_type = 'cells';
        tuple.calcium_trace = single(mean(cellTraces,2));
        tuple.annulus_trace = single(mean(annulusTraces,2));
        if redChannel;tuple.red_trace = single(mean(redTraces,2));end
        insert( obj, tuple );
        
        % insert site entry with masknum=-3
        tuple=key;
        tuple.masknum=-3;
        tuple.mask_type = 'site';
        tuple.calcium_trace = single(mean([cellTraces neuropilTrace],2));
        if redChannel;tuple.red_trace = single(mean([redTraces neuropilrTrace],2));end
        insert( obj, tuple );
        
        % insert cells into table, one at a time
        for iCell = 1:length(x);
            tuple=key;
            tuple.masknum = iCell;
            tuple.calcium_trace = single(cellTraces(:,iCell));
            tuple.annulus_trace = single(annulusTraces(:,iCell));
            if redChannel;tuple.red_trace     = single(redTraces(:,iCell));end
            if ~isNeuron(iCell) && ~isempty(strfind(dyes,'SR'))
                tuple.mask_type = 'astrocyte';
            else
                tuple.mask_type = 'neuron';
            end
            insert( obj, tuple );
            tuple=key;
            tuple.masknum = iCell;
            tuple.green_contrast= greenContrast(iCell);
            if redChannel;tuple.red_contrast  = redContrast(iCell);end
            tuple.img_x = x(iCell);
            tuple.img_y = y(iCell);
            tuple.cell_radius = radii(iCell);
            tuple.sharpness = sharpness(iCell);
            tuple.fit_qual = correlation(iCell);  % correlation between the fitted blob and the pixel values
            makeTuples(MaskCells,tuple)
        end
        
        % insert patched cell with masknum=0
        if  strcmp(fetch1(Scans(key),'aim'),'patching')
            
            % get ephys data
            elChan  = readCh4(tpReader(Scans(key)));
            eData = elChan(:);
            tuple = key;
            tuple.masknum = 0;
            tuple.mask_type = 'ephys';
            tuple.calcium_trace = [];
            tuple.red_trace = [];
            tuple.ephys_trace = eData;
            insert( obj, tuple );
        end
    end
end

%% MPScan
if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
    % get site data
    dyes = fetch1(Experiments(key),'dyes');
    [greenImg,redImg] = fetch1( Movies(key),'affine_green','affine_red');
    tpr = tpReader(Scans(key));
    if isempty(greenImg) || isempty(redImg)
        disp('Aligned images are missing, computing on non aligned data');
        [greenImg,redImg] = fetch1( Movies(key),'raw_green','raw_red');
        imChan = getImagingChannel(tpr,1);
        imChan2 = getImagingChannel(tpr,2);
    else
        % extract traces
        fileNameCore = getLocalPath([fetch1(Experiments(key),'directory') '/' fetch1(Scans(key),'file_name')]);
        if exist(getLocalPath([fileNameCore '_fineAlignment.mat']),'file')
            channel1 = getImagingChannel(tpr,1);
            channel2 = getImagingChannel(tpr,2);
            % get the aligning info
            load(getLocalPath([fileNameCore '_fineAlignment.mat']));
            % apply alignment
            degrees = [mc.warp_degree mc.warp_degree 0];
            imChan = zeros(size(channel1),'uint16');
            for i=1:size(channel1,3)
                imChan(:,:,i) = uint16(tpMethods.YWarp.apply(double(channel1(:,:,i)), mc.warp_polynom(i,:), degrees));
            end
            imChan2 = zeros(size(channel2),'uint16');
            for i=1:size(channel2,3)
                imChan2(:,:,i) = uint16(tpMethods.YWarp.apply(double(channel2(:,:,i)), mc.warp_polynom(i,:), degrees));
            end
        elseif isAligned(tpr)
            imChan = getAlignedImagingChannel(tpr,1);
            imChan2 = getAlignedImagingChannel(tpr,2);
        elseif ~isAligned(tpr)
            % get the aligning info
            fileName = [fetch1(Scans(key),'file_name') '_alignment.mat'];
            path = fetch1(Experiments(key),'directory');
            alignData = load(getLocalPath([path '/' fileName]));
            
            % apply alignment
            stack = getData(getImagingChannel(tpr,1));
            stack = transformStack(stack, alignData.alignData.transMatrices);
            validRect = alignData.alignData.validRect;
            imChan = stack(validRect(3):validRect(4), validRect(1):validRect(2), :);
            
            stack = getData(getImagingChannel(tpr,2));
            stack = transformStack(stack, alignData.alignData.transMatrices);
            validRect = alignData.alignData.validRect;
            imChan2 = stack(validRect(3):validRect(4), validRect(1):validRect(2), :);
        end
    end
    if nargin>2 % manualy insert data
        [cellTraces,~,annulusTraces] = getTraces( cellFinder, imChan );
        redTraces = getTraces( cellFinder, imChan2 );
        
        tuple=key;
        tuple.calcium_trace = single(cellTraces);
        tuple.annulus_trace = single(annulusTraces);
        tuple.red_trace     = single(redTraces);
        if nargin>3
            tuple.mask_type = type;
        else
            tuple.mask_type = 'neuron';
        end
        insert( obj, tuple );
        
        tuple=key;
        [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder ); %#ok<ASGLU>
        greenContrast = getCellContrast( cellFinder, greenImg );
        redContrast   = getCellContrast( cellFinder,   redImg );
        tuple.green_contrast= greenContrast;
        tuple.red_contrast  = redContrast;
        tuple.img_x = x;
        tuple.img_y = y;
        tuple.cell_radius = radii;
        tuple.sharpness = sharpness;
        tuple.fit_qual = correlation;  % correlation between the fitted blob and the pixel values
        makeTuples(MaskCells,tuple)
        
    else
        % fit the green image as a sum of radially symmetric blobs
        % with radii proportional to the total magnification
        % find cells
        [lens,mag] = fetch1(Scans(key),'lens','mag');
        pixelPitch = 11000/512/mag/lens;
        
        cellFinder = CellFinder( greenImg, pixelPitch, 'minRadius', 3.0, 'minDistance', 7.0, 'minContrast', 0.05 );
        fig = figure(1000);
        plot(cellFinder)
        print(fig,'-dpng', [fetch1(Experiments(key),'directory') '\'...
            fetch1(Scans(key),'file_name') '_fcells'])
        
        
        [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder ); %#ok<ASGLU>
        greenContrast = getCellContrast( cellFinder, greenImg );
        redContrast   = getCellContrast( cellFinder,   redImg );
        isNeuron  =  2.5 > (redContrast-trimmean(redContrast,50))/iqr(redContrast)/0.7413;  % rule of thumb for neurons
        
        [cellTraces,neuropilTrace,annulusTraces] = getTraces( cellFinder, imChan );
        
        [redTraces,neuropilrTrace] = getTraces( cellFinder, imChan2 );
        
        % insert neuropil entry with masknum=-1
        tuple=key;
        tuple.masknum=-1;
        tuple.mask_type = 'neuropil';
        tuple.calcium_trace = single(neuropilTrace);
        tuple.red_trace = neuropilrTrace;
        insert( obj, tuple );
        
        % insert cells entry with masknum=-2
        tuple=key;
        tuple.masknum=-2;
        tuple.mask_type = 'cells';
        tuple.calcium_trace = single(mean(cellTraces,2));
        tuple.annulus_trace = single(mean(annulusTraces,2));
        tuple.red_trace = single(mean(redTraces,2));
        insert( obj, tuple );
        
        % insert site entry with masknum=-3
        tuple=key;
        tuple.masknum=-3;
        tuple.mask_type = 'site';
        tuple.calcium_trace = single(mean([cellTraces neuropilTrace],2));
        tuple.red_trace = single(mean([redTraces neuropilrTrace],2));
        insert( obj, tuple );
        
        % insert cells into table, one at a time
        for iCell = 1:length(x);
            tuple=key;
            tuple.masknum = iCell;
            tuple.calcium_trace = single(cellTraces(:,iCell));
            tuple.annulus_trace = single(annulusTraces(:,iCell));
            tuple.red_trace     = single(redTraces(:,iCell));
            if ~isNeuron(iCell) && ~isempty(strfind(dyes,'SR'))
                tuple.mask_type = 'astrocyte';
            else
                tuple.mask_type = 'neuron';
            end
            insert( obj, tuple );
            tuple=key;
            tuple.masknum = iCell;
            tuple.green_contrast= greenContrast(iCell);
            tuple.red_contrast  = redContrast(iCell);
            tuple.img_x = x(iCell);
            tuple.img_y = y(iCell);
            tuple.cell_radius = radii(iCell);
            tuple.sharpness = sharpness(iCell);
            tuple.fit_qual = correlation(iCell);  % correlation between the fitted blob and the pixel values
            makeTuples(MaskCells,tuple)
        end
        
        % insert patched cell with masknum=0
        if  strcmp(fetch1(Scans(key),'aim'),'patching')
            
            % get ephys data
            elChan  = getElectrodeChannel(tpr,2);
            eData = elChan(:);
            tuple = key;
            tuple.masknum = 0;
            tuple.mask_type = 'ephys';
            tuple.calcium_trace = [];
            tuple.red_trace = [];
            tuple.ephys_trace = eData;
            insert( obj, tuple );
        end
    end
end

%% AOD Scan
if strcmp(fetch1(Scans(key),'scan_prog'),'AOD')
    
    [path,name] = fetch1( Experiments*Scans(key), 'directory','file_name' );
    fps = fetch1(Movies(key),'fps');
    
    if ~isempty(Scans(key,'exp_date > "2012-05-01"'))
        dirs = dir(['M:\Mouse\' name(1:10) '*']);
        names = vertcat(dirs.name);
        timediff = str2num(names(:,[12 13 15 16 18 19]))- str2num( name([12 13 15 16 18 19])); %#ok<ST2NM>
        timeind = find(timediff<0);
        [~,i] = max(timediff(timeind));
        itime = timeind(i);
        filename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' name];
        [traces, ~, coordinates, traces2] = loadAODTraces(filename,fps);
    else
        filename = getLocalPath([path '/' name '.h5']);
        [traces, ~, coordinates, traces2] = loadTraces(filename,fps);
    end
    
    % detect & remove duplicates
    thr = 5; % threshold in microns
    mtraces = mean(traces); % get mean of traces
    d = squareform(pdist(coordinates)); % get pairwise distances
    d(logical(tril(ones(size(d,1),size(d,2)),0)))=thr; % select upper triangle
    [indx,indy] = ind2sub(size(d),find(d<thr));ind=[indx indy];% get indexes
    [~,sind]=sort(mtraces(ind),2,'ascend'); % select low firing ones
    traces(:,ind(sind==1))      = []; % remove extra cells
    traces2(:,ind(sind==1))     = []; % remove extra cells
    coordinates(ind(sind==1),:) = []; % remove extra cells
    
    for itrace = 1:size(traces,2)
        tuple=key;
        tuple.masknum = itrace;
        tuple.mask_type = 'neuron';
        tuple.calcium_trace = single(traces(:,itrace));
        tuple.red_trace = single(traces2(:,itrace));
        insert( obj, tuple );
        
        tuple=key;
        tuple.masknum = itrace;
        tuple.green_contrast = 0;
        tuple.red_contrast  = 0;
        tuple.img_x = coordinates(itrace,1);
        tuple.img_y = coordinates(itrace,2);
        tuple.img_z = coordinates(itrace,3);
        tuple.cell_radius = 0;
        tuple.sharpness = 0;
        tuple.fit_qual = 0;  % correlation between the fitted blob and the pixel values
        makeTuples(MaskCells,tuple)
    end
end

%% Unirec Scan
if strcmp(fetch1(Scans(key),'scan_prog'),'Unirec')
    
    tuple=key;
    tuple.masknum = 1;
    tuple.mask_type = 'multiunit';
    tuple.calcium_trace = [];
    tuple.red_trace = [];
    insert( obj, tuple );
    
    tuple=key;
    tuple.masknum = 1;
    tuple.green_contrast = 0;
    tuple.red_contrast  = 0;
    tuple.img_x = 0;
    tuple.img_y = 0;
    tuple.img_z = 0;
    tuple.cell_radius = 0;
    tuple.sharpness = 0;
    tuple.fit_qual = 0;  % correlation between the fitted blob and the pixel values
    makeTuples(MaskCells,tuple)
end