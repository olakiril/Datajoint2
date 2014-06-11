function fh = plot(keys)


keys = fetch(keys);
fh = figure;
for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
    if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        clf
        
        % Load movie information
        tpr = tpReader(Scans(key));
        % detect cells
        imChan = getImagingChannel(tpr,1);
        imChan2 = getImagingChannel(tpr,2);
        
        [lens,mag] = fetch1(Scans(key),'lens','mag');
        pixelPitch = 11000/512/mag/lens;
        
        steps = fetch1(Stacks(key),'steps');
        for istack = 1:steps
            
            if strcmp(fetch1(Experiments(key),'dyes'),'SR')
                greenImg = double(abs(imChan2(:,:,istack) - max(max(imChan2(:,:,istack)))));
                % enchance image
                greenImg(greenImg < prctile(greenImg(:),1)) = prctile(greenImg(:),1);
                greenImg(greenImg > prctile(greenImg(:),99)) = prctile(greenImg(:),99);
                
            else
                greenImg = double(imChan(:,:,istack));
            end
            
            [x y radii sharpness contrast ] = ...
                fetchn(MaskCells(key,['masknum>0 and img_z = ' num2str(istack)]),...
                'img_x','img_y','cell_radius','sharpness','green_contrast');
            
            if isempty(x)
                
                cellFinder = CellFinder( greenImg, pixelPitch,'minRadius', 3.0, 'minDistance', 7.0, ...
                    'minContrast',0.001,'minSharpness',3,'minCorr',0.5);
                
                [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder );
                
                im = normalize(greenImg);
                im = im(1:floor(numel(greenImg)/16)*16);
                im = reshape(im,16,[]);
                s = mean(std(im));
                a1 =       5.114 ;
                b1 =      0.1854 ;
                c1 =     0.06138 ;
                coef =   a1*exp(-((s-b1)/c1)^2);
                indx = (correlation>median(correlation)-coef*std(correlation)) &...
                    (contrast>median(contrast)-std(contrast)/4);
                x = x(indx);
                y = y(indx);
                radii = radii(indx);
                contrast = contrast(indx);
                sharpness = sharpness(indx);
            end
            
            cellFinder = CellFinder( greenImg, pixelPitch,...
                'x', x, 'y', y, 'radii', radii,'contrast',contrast,...
                'sharpness',sharpness,'cleanImg',greenImg);
            
            plot(cellFinder,'color',[0.3 0.3 0.3],'text',0,'plotsmall',1)
            pause
        end
    end
end