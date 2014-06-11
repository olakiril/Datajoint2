function obj = compute( obj, fieldname, key )

fprintf( 'Computing %s/%s...\n', class( obj ), fieldname );

switch fieldname
    case 'job'
        obj = update( obj, 'job', 0 );
        temp = take( CircTraces(key), 'skewness' );
        obj = update( obj, 'job', 1 );

    case {'nframes','fps','total_size','first_valid_frame','last_valid_frame'}
        tpr = tpReader(obj);
        channel = getImagingChannel(tpr,1);
        nframes = size( channel, 3 );
        totalSize = prod(size(channel));
        fps = getFrameRate( tpr );
        times = readTimestamps( tpr );
        if size(times,2)==1, times=times', end;
        minDist = 2.5*median(diff(times));  % frames must not be separated by more than 2 frame periods
        fc = conseqCount( times, minDist );  % forward count
        bc = fliplr( conseqCount( fliplr(times), minDist ) );
        validFrames = find( (fc+bc)==max(fc+bc) );

        obj = update( obj, 'nframes', nframes, 'fps', fps, 'total_size', totalSize...
            , 'first_valid_frame', validFrames(1), 'last_valid_frame', validFrames(end) );

    case {'raw_green','raw_red'}
        tpr = tpReader(obj);
        channel = getImagingChannel(tpr,1+strcmp(fieldname,'raw_red'));
        block = 128;   %  compromise between speed and memory performance
        nFrames = size(channel,3);
        img = 0;
        for iFrame=1:block:nFrames
            img = img + sum(double(channel(:,:,iFrame:min(iFrame+block-1,nFrames))),3);
        end
        obj = update( obj, fieldname, single(img/nFrames) );

    case {'affine_green','affine_red'}
        tpr = tpReader(obj);
        channel = getAlignedImagingChannel(tpr,1+strcmp(fieldname,'affine_red'));
        block = 128;   %  compromise between speed and memory performance
        nFrames = size(channel,3);
        img = 0;
        for iFrame=1:block:nFrames
            img = img + sum(double(channel(:,:,iFrame:min(iFrame+block-1,nFrames))),3);
        end
        obj = update( obj, fieldname, single(img/nFrames) );

    case {'green','red'}
        tpr = tpReader(obj);
        channel = getImagingChannel(tpr,1+strcmp(fieldname,'red'));
        [rasterCorrection,motionCorrection] = take( obj, 'raster_correction', 'motion_correction' );
        block = 128;   %  compromise between speed and memory performance
        nFrames = size(channel,3);
        img = 0;
        for iFrame=1:block:nFrames
            m = double(channel(:,:,iFrame:min(iFrame+block-1,nFrames)));
            for jFrame=1:size(m,3)
                m(:,:,jFrame) = fixInterlace( m(:,:,jFrame), 'initialSolution', rasterCorrection(iFrame+jFrame-1,:), 'maxIter', 0);
                m(:,:,jFrame) = alignFrame( m(:,:,jFrame), motionCorrection(iFrame+jFrame-1,:) );
            end
            img = img + sum(m,3);
        end
        obj = update( obj, fieldname, single(img/nFrames) );

    case 'raster_correction'
        % compute raster correction based on green images only
        chan = getImagingChannel(tpReader(obj), 1);
        sz = size(chan);

        % perform optimization iterations for every 'step' frames
        step = 60;  %must be even
        maxIter = 10;  % number of optimization iterations for the seed solution
        oddCorrection = [];
        evenCorrection = [];
        if sz(3) > 0
            for iFrame = 1:step:sz(3)
                if iFrame>sz(3)-2*step+1
                    step = sz(3)-iFrame+1;
                end
                fprintf('[%4d/%4d]\n',iFrame,sz(3));
                oddFrames  = iFrame+0:2:iFrame+step-1;
                evenFrames = iFrame+1:2:iFrame+step-1;

                img = mean(anscombe(chan(:,:,oddFrames)),3);
                [img,oddCorrection] = fixInterlace( img, 'initialSolution',oddCorrection, 'maxIter',maxIter);
                for k=oddFrames
                    rasterCorrection(k,:) = oddCorrection;
                end

                img = mean(anscombe(chan(:,:,evenFrames)),3);
                [img,evenCorrection] = fixInterlace( img,'initialSolution',evenCorrection, 'maxIter',maxIter);
                for k=evenFrames
                    rasterCorrection(k,:) = evenCorrection;
                end
                maxIter = 2;  % update solution using a smaller number of iterations (to save time and avoid overfitting)

            end
            obj = update( obj, 'raster_correction', rasterCorrection );
        end



    case 'motion_correction'
        % compute motion correction based on green images
        chan = getImagingChannel(tpReader(obj),1);
        sz = size( chan );
        m = anscombe(chan(:,:,1:min(sz(3),30)));
        rasterCorrection = take( obj, 'raster_correction' );
        for iFrame=1:size(m,3)
            m(:,:,iFrame)=fixInterlace(m(:,:,iFrame),'initialSolution',rasterCorrection(iFrame,:),'maxIter',0);
        end
        refFrame = median( m, 3 );

        block = 50;
        for iFrame = 1:sz(3)
            if mod(iFrame,block)==1
                m = anscombe(chan(:,:,iFrame:min(sz(3),iFrame+block-1)));
                fprintf('[%4d/%4d]\n',iFrame-1,sz(3));
            end
            img = m(:,:,mod(iFrame-1,block)+1);
            img = fixInterlace( img, 'initialSolution', rasterCorrection(iFrame,:), 'maxIter', 0 );
            p = alignFrame( img, refFrame );
            img = alignFrame( img, p );
            refFrame = refFrame + (img-refFrame)/50;  % refine refFrame
            motionCorrection(iFrame,:)=p;
            imagesc( img, [0.05 0.65] ); axis image; disp(p);
            drawnow;
        end
        obj = update( obj, 'motion_correction', motionCorrection );


    otherwise
        warning( '%s/compute does not know how to compute %s', class( obj ), fieldname );
end



