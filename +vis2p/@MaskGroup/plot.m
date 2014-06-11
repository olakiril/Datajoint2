function fh = plot(keys)

keys = fetch(keys);
fh = figure;
% set(gcf,'position',[ 11   260   889   350]);
% set(gcf,'position',[ 2022         249        1463         719])
clf
for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
    if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        %% MPScan
        clf
        
        greenImg = fetch1( Movies(key),'affine_green');
        if isempty(greenImg)
            greenImg = fetch1( Movies(key),'raw_green');
        end
        
        [lens,mag] = fetch1(Scans(key),'lens','mag');
        pixelPitch = 11000/512/mag/lens;
        
        [x, y, radii, sharpness, contrast ] = ...
            fetchn(MaskCells(key,'masknum>0'),...
            'img_x','img_y','cell_radius','sharpness','green_contrast');
        
        obj = CellFinder( greenImg, pixelPitch,...
            'x', x, 'y', y, 'radii', radii,'contrast',contrast,'sharpness',sharpness,'cleanImg',greenImg/max(greenImg(:)));
        
        types = fetchn(MaskTraces(key,'masknum>0'),'mask_type');
        colors = zeros(length(types),3);
        colors(cellfun(@(x) strcmp(x,'neuron'),types),3) = 1;
        colors(cellfun(@(x) ~strcmp(x,'neuron'),types),1) = 1;
        %         plot(cellFinder)
        [gr, rd] = fetch1(Movies(key),'affine_green','affine_red');
        subplot(221)
        imshow(normalize(gr))
        hold on;
        for i=1:length(x)
            rectangle( 'Position', [[x(i),y(i)]-radii(i), 2*radii(i)*[1 1]], 'Curvature', [1 1],'EdgeColor',colors(i,:));
            text(x(i),y(i),int2str(i), 'Color',[1 0.2 0]);
        end
        hold off
        subplot(223)
        imshow(normalize(gr))
        
        subplot(222)
%         imshow(normalize(normalize(rd)-normalize(gr)))
        imshow(normalize(normalize(rd)))
         hold on;
        for i=1:length(x)
            rectangle( 'Position', [[x(i),y(i)]-radii(i), 2*radii(i)*[1 1]], 'Curvature', [1 1],'EdgeColor',colors(i,:));
            text(x(i),y(i),int2str(i), 'Color',[1 0.2 0]);
        end
        hold off
        subplot(224)
        
        imshow(normalize(normalize(rd)-normalize(gr)))
        
        if ikey ~= 1
            pause
        end
        
    elseif strcmp(fetch1(Scans(key),'scan_prog'),'ScanImage')
        %% ScanImage
        clf
        
        greenImg = fetch1( Movies(key),'affine_green');
        if isempty(greenImg)
            greenImg = fetch1( Movies(key),'raw_green');
        end
        
        [lens,mag] = fetch1(Scans(key),'lens','mag');
        pitch = 11000/mag/lens;
        pixelPitch = pitch/size(greenImg,2);
          
        [x, y, radii, sharpness, contrast ] = ...
            fetchn(MaskCells(key,'masknum>0'),...
            'img_x','img_y','cell_radius','sharpness','green_contrast');
        
        obj = CellFinder( greenImg, pixelPitch,...
            'x', x, 'y', y, 'radii', radii,'contrast',contrast,'sharpness',sharpness,'cleanImg',greenImg/max(greenImg(:)));
        
        types = fetchn(MaskTraces(key,'masknum>0'),'mask_type');
        colors = zeros(length(types),3);
        colors(cellfun(@(x) strcmp(x,'neuron'),types),3) = 1;
        colors(cellfun(@(x) ~strcmp(x,'neuron'),types),1) = 1;
                plot(obj)
        [gr, rd] = fetch1(Movies(key),'affine_green','affine_red');
        
        if isempty(rd)
            [rastcor, motcor] = fetch1(Movies(key),'raster_correction','motion_correction');
            rd = mean(getFrames(Movies(key), 2,1:33,rastcor, motcor),3);
        end
        
        subplot(221)
        imshow(normalize(gr))
               hold on;
        for i=1:length(x)
            rectangle( 'Position', [[x(i),y(i)]-radii(i), 2*radii(i)*[1 1]], 'Curvature', [1 1],'EdgeColor',colors(i,:));
            text(x(i),y(i),int2str(i), 'Color',[1 0.2 0]);
        end
        hold off
        subplot(223)
        imshow(normalize(gr))
        im = normalize(rd);
      clim = [min(min(im(70:end,:))) max(max(im(70:end,:)))];
        subplot(222)
        imshow(im,clim)
         hold on;
        for i=1:length(x)
            rectangle( 'Position', [[x(i),y(i)]-radii(i), 2*radii(i)*[1 1]], 'Curvature', [1 1],'EdgeColor',colors(i,:));
            text(x(i),y(i),int2str(i), 'Color',[1 0.2 0]);
        end
        hold off
        subplot(224)
        imshow(im,clim)
        
        if ikey ~= 1
            pause
        end
    end
end