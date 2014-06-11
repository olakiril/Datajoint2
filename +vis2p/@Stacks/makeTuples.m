function makeTuples( obj, key )

%% MPScan
if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
    
    % Load movie information
    tpr = tpReader(Scans(key));
    
    % insert stack properties
    properties = getProperties(tpr);
    tuple = key;
    tuple.interval = properties.intervalZ;
    tuple.steps = properties.sectionCount;
    tuple.averages = properties.averagingCount;
    tuple.xsize = properties.frameSize(1);
    tuple.ysize = properties.frameSize(2);
    
%     insert( obj, tuple );
    
    % detect cells
    imChan = getImagingChannel(tpr,1);
    imChan2 = getImagingChannel(tpr,2);
    
    [lens,mag] = fetch1(Scans(key),'lens','mag');
    pixelPitch = 11000/512/mag/lens;
    
    maskend = 0;
    for iStep = 1:size(imChan,3);
        try
            if strcmp(fetch1(Experiments(key),'dyes'),'SR')
                cimage = abs(imChan2(:,:,iStep) - max(max(imChan2(:,:,iStep))));
                cimage(cimage < prctile(cimage(:),1)) = prctile(cimage(:),1);
                cimage(cimage > prctile(cimage(:),99)) = prctile(cimage(:),99);
                
            else
                cimage = imChan(:,:,iStep);
            end
            
            cellFinder = CellFinder( cimage, pixelPitch,'minRadius', 3.0, 'minDistance', 7.0, ...
                'minContrast',0.001,'minSharpness',3,'minCorr',0.5);
            
            [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder ); 
            
            redContrast   = getCellContrast( cellFinder,   imChan2(:,:,iStep) );
            greenContrast   = getCellContrast( cellFinder,   imChan(:,:,iStep) );
            
            % select nicely fitted cells
            im = normalize(cimage);
            im = im(1:floor(numel(cimage)/16)*16);
            im = reshape(im,16,[]);
            
            s = mean(std(double(im)));
            a1 =       5.114 ;
            b1 =      0.1854 ;
            c1 =     0.06138 ;
            coef =   a1*exp(-((s-b1)/c1)^2);
            indx = (correlation>median(correlation)-coef*std(correlation)) & (contrast>median(contrast)-std(contrast)/4);
            display([num2str(sum(indx)) '/' num2str(length(x)) ' passed increased threshold'])
            x = x(indx); y = y(indx); radii = radii(indx);
            sharpness = sharpness(indx);correlation = correlation(indx);
            redContrast = redContrast(indx); greenContrast = greenContrast(indx);
            
            % insert cells into table, one at a time
            for iCell = 1:length(x);
                if iCell ==1
                    display(['Inserting ' num2str(length(x)) ' Cells'])
                end
                tuple=key;
                tuple.masknum = iCell + maskend;
                tuple.green_contrast= greenContrast(iCell);
                tuple.red_contrast  = redContrast(iCell);
                
                if isnan(tuple.green_contrast) || isinf(tuple.green_contrast)
                    tuple.green_contrast = 0;
                end
                if isnan(tuple.red_contrast) || isinf(tuple.red_contrast)
                    tuple.red_contrast = 0;
                end
                
                tuple.img_x = x(iCell);
                tuple.img_y = y(iCell);
                tuple.img_z = iStep;
                tuple.cell_radius = radii(iCell);
                tuple.sharpness = sharpness(iCell);
                tuple.fit_qual = correlation(iCell);  % correlation between the fitted blob and the pixel values
                makeTuples(MaskCells,tuple)
            end
            
            if ~isempty(x)
                maskend = tuple.masknum;
            end
        catch
            display(['Skipping step: ' num2str(iStep)])
            continue
        end
    end
end
