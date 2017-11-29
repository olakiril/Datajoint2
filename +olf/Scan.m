%{
olf.Scan (imported) #
-> olf.Session
---
fps                         : double                        # frames per second
raw_ch1                     : mediumblob                    # mean ch1 image before any corrections
raw_ch2                     : mediumblob                    # mean ch2 image before any corrections
nframes                     : smallint unsigned             # the number of frames
xsize                       : smallint                      # scan x pixel size
ysize                       : smallint                      # scan y pixel size
zsize                       : smallint                      # scan z slice size
fov=null                    : int                           # Field of View in microns
x=null                      : int                           # m) (um) objective manipulator x position
y=null                      : int                           # m) (um) objective manipulator y position
z=null                      : int                           # m) (um) objective manipulator z position
trials                      : mediumblob                    # decoded trial time
%}


classdef Scan < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = (olf.Session)
    end
    
    methods (Access=protected)
        function makeTuples(obj,key)
            
            [file, path, mag] = fetch1(olf.Session(key),'file_name','path','lens');
            filetype = getLocalPath(fullfile(path,sprintf('%s%s',file,'*')));
            display(['Reading file: ' filetype])
            file_dirs = dir(filetype);
            tuple=key;
            dur = fetch1(olf.Session(key),'pulse_duration');
            
            if  file_dirs(1).bytes> 2*10^9 % hack to read large files
                tr = readtiffinfo(obj,getLocalPath(fullfile(path,file_dirs(1).name)));
                pos = tr.scanimage.SI.hMotors.motorPosition;
                
                tuple.xsize = tr.scanimage.SI.hRoiManager.pixelsPerLine;
                tuple.ysize = tr.scanimage.SI.hRoiManager.linesPerFrame;
                tuple.zsize = tr.scanimage.SI.hStackManager.numSlices; % fix this
                if tuple.zsize>1
                    tuple.fps = tr.scanimage.SI.hRoiManager.scanVolumeRate;
                else
                    tuple.fps = tr.scanimage.SI.hRoiManager.scanFrameRate;
                end
                zoom = tr.scanimage.SI.hRoiManager.scanZoomFactor;
                channels = length(tr.scanimage.SI.hChannels.channelSave);
                
                % get Data
                data_ch1 = int16([]);
                data_ch2 = int16([]);
                data_phd = int16([]);
                if tuple.zsize>1
                    for ifile = 1:length(file_dirs);
                        data = read_patterned_tifdata(getLocalPath(fullfile(path,file_dirs(ifile).name)));
                        chunkSize = size(data,3)/(channels*tuple.zsize);
                        idx = 1+(ifile-1)*chunkSize:(ifile)*chunkSize;
                        for islice = 1:tuple.zsize
                            data_ch1(:,:,idx,islice) = ...
                                data(:,:,1 + channels*(islice-1):channels*tuple.zsize:end);
                            data_phd(:,:,idx,islice) = ...
                                data(:,:,3 + channels*(islice-1):channels*tuple.zsize:end);
                        end
                        clear data
                    end
                    data_phd = permute(data_phd,[1 2 4 3]);
                else
                    for ifile = 1:length(file_dirs);
                        data = read_patterned_tifdata(getLocalPath(fullfile(path,file_dirs(ifile).name)));
                        data_ch1(:,:,end+1:end+size(data,3)/3) = data(:,:,1:channels*tuple.zsize:end);
                        data_ch2(:,:,end+1:end+size(data,3)/3) = data(:,:,2:channels*tuple.zsize:end);
                        data_phd(:,:,end+1:end+size(data,3)/3)= data(:,:,3:channels:end);
                        clear data
                    end
                end
            else
                reader = ne7.scanimage.Reader5(filetype);
                pos = reader.header.hMotors_motorPosition(1:3);
                sz = size(reader);
                tuple.fps = reader.fps;
                tuple.xsize = sz(1);
                tuple.ysize = sz(2);
                tuple.zsize = sz(4);
                tuple.nframes = sz(5);
                zoom = reader.zoom;
                
                % get Data
                data_ch1 = permute(squeeze(reader(:,:,1,:,:)),[1 2 4 3]);
                data_phd = reader(:,:,end,:,:);
                
            end
            
            % filter photodiode
            tuple.trials = findTrials(obj,data_phd,tuple.fps,dur,tuple.zsize);
            clear data_phd
            
            % fix baseline
            data_ch1 = data_ch1+abs(min(data_ch1(:)));
            
            % fix raster & motion artifacts
            for islice = 1:size(data_ch1,4)
                if true % motion correction
                    dat = squeeze(data_ch1(:,:,:,islice));
                    offsets = tpMethods.MotionCorrection.fit(single(dat));
                    offsets = bsxfun(@minus, offsets, median(offsets));
                    data_ch1(:,:,:,islice) = tpMethods.MotionCorrection.apply(dat, offsets);
                end
            end
            
            tuple.nframes = size(data_ch1,3);
            tuple.x = pos(1); tuple.y = pos(2);tuple.z = pos(3);
            tuple.fov = 21000/mag/zoom;
            tuple.raw_ch1 = squeeze(mean(data_ch1,3));
            tuple.raw_ch2 = zeros(size(tuple.raw_ch1));
            
            % insert Scan
            insert( obj, tuple );
            
            % extract Cell Traces
            pixelPitch = tuple.fov/size(tuple.raw_ch1,2);
            
            if false % auto segmentation
                cellFinder = CellFinder( tuple.raw_ch1(:,:,1), pixelPitch, 'minRadius', 3.0, 'minDistance', 7.0, 'minContrast', 0.5 );
                figure
                plot(cellFinder)
                [x,y,radii,contrast,sharpness,correlation] = getCells( cellFinder );
                [cellTraces,neuropilTrace] = getTraces( cellFinder, data_ch1 );
                [redTraces,neuropilrTrace] = getTraces( cellFinder, data_ch2 );
                masks = getMasks(cellFinder);
            else % manual segmentation
                data_ch1 = permute(data_ch1,[3 1 2 4]);
                
                % initialize variables
                neuropilTrace = zeros(size(data_ch1,1),1);
                mask_idx = 0;
                cellTraces = [];
                masks = [];radii = [];x = [];y = [];z = [];
                for islice = 1:tuple.zsize
                    data = data_ch1(:,:,:,islice);
                    logimg = log(max(max(tuple.raw_ch1(:)/1000),tuple.raw_ch1(:,:,islice)));
                    k = hamming(2*round(2.0*8/pixelPitch)+1);  k=k/sum(k);  %  subtract the background
                    sharpImg = logimg-imfilter(imfilter(logimg,k,'symmetric'),k','symmetric');
                    fig = figure;
                    imagesc(sharpImg);axis image;colormap gray
                    imCells = drawCells;close(fig)
                    mask_labels = bwlabel(imCells,4);
                    unique_masks = unique(mask_labels(:));
                    
                    neuropilTrace = neuropilTrace + mean(data(:,mask_labels == 0),2);
                    
                    for imask = 1:length(unique_masks)-1
                        mask = mask_labels == imask;
                        mask_idx = mask_idx+1;
                        stats = regionprops('table',mask,'Centroid',...
                            'MajorAxisLength','MinorAxisLength');
                        x(mask_idx) = stats.Centroid(1);
                        y(mask_idx) = stats.Centroid(2);
                        z(mask_idx) = islice;
                        diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
                        radii(mask_idx) = diameters/2;
                        cellTraces(:,mask_idx) = mean(data(:,mask),2);
                        masks{mask_idx} = mask;
                    end
                end
                neuropilTrace = neuropilTrace/tuple.zsize;
                redTraces = zeros(size(cellTraces));
                neuropilrTrace = zeros(size(neuropilTrace));
            end
            
            % insert mask information
            for imask = 0:1:length(masks)
                if imask==0;
                    mask_type = 'neuropil';
                    trace_ch1 = single(neuropilTrace);
                    trace_ch2 = single(neuropilrTrace);
                    xpos = 0;
                    ypos = 0;
                    zpos = 0;
                    rad = 0;
                    mask = ones(size(tuple.raw_ch1));
                    for imask0 = 1:length(masks)
                        mask(masks{imask0}) = 0;
                    end
                else
                    mask_type='neuron';
                    trace_ch1 = single(cellTraces(:,imask));
                    trace_ch2 = single(redTraces(:,imask));
                    xpos = x(imask);
                    ypos = y(imask);
                    zpos = z(imask);
                    rad = radii(imask);
                    mask = zeros(size(tuple.raw_ch1));
                    mask( masks{imask}) = 1;
                end
                makeTuples(olf.MaskTraces,key, imask, mask, mask_type, trace_ch1, trace_ch2,xpos,ypos,zpos,rad)
            end
        end
    end
    
    methods
        function self = Scan(varargin)
            self.restrict(varargin{:})
        end
        
        function tif_info = readtiffinfo(obj, filename)
            tr = Tiff(filename);
            id = getTag(tr,'ImageDescription');
            d = find(id==10);
            for i = 2:length(d)-1;
                try
                    eval(sprintf('%s%s',id(d(i-1)+1:d(i)-1),';'))
                end
                
            end
            clear tr; clear id; clear d;clear filename
            w = whos;
            for a = 1:length(w)
                tif_info.(w(a).name) = eval(w(a).name);
            end
            tif_info = rmfield(tif_info,'obj');
        end
        
        function trials = findTrials(obj,data_phd,fps,dur,nslices)
            
            phd = reshape(mean(data_phd,2),[],1);
            phd = (mean(phd)+min(phd))/2-phd > 0;
            
            flipdur = nslices*dur*fps*size(data_phd,1)/33/1000; % in line time
            all_times = find(diff(phd));
            start_times = find(diff(phd)>0);
            prod = 2 .^ (0 : 15 )';
            istart = 1;
            while istart <= length(all_times)
                lookuptimes = all_times(istart)+flipdur/2+flipdur*(1:2:32);
                b = phd(round(lookuptimes))';
                trial = b*prod;
                next = find(start_times>lookuptimes(end)+flipdur*3/2,1,'first');
                stimtimes = 1:size(phd);
                stim = zeros(size(stimtimes));
                stim(stimtimes > all_times(istart) & stimtimes < all_times(istart)+33*flipdur) = trial;
                previous_trial = trial;
                while start_times(next) < all_times(end)
                    
                    nextlookuptimes = all_times(start_times(next)==all_times)+flipdur/2+flipdur*(1:2:32);
                    b = phd(round(nextlookuptimes))';
                    next_trial = b*prod;
                    stim(stimtimes > all_times(start_times(next)==all_times) ...
                        & stimtimes < all_times(start_times(next)==all_times)+33*flipdur) = next_trial;
                    
                    next = find(start_times>nextlookuptimes(end)+flipdur*2,1,'first');
                    if next_trial ~= previous_trial + 1
                        break
                    elseif start_times(next) < all_times(end)
                        istart = length(all_times)+1;
                    end
                    previous_trial = next_trial;
                end
                istart = istart+ 1;
            end
            trials = stim(1:size(data_phd,1):end);
        end
        
        function plotMask(obj)
            [fov,im] = fetch1(obj,'fov','raw_ch1');
            
            figure
            imagesc(im);hold on
            colormap gray
            axis image
            axis off
            plot([size(im,1)*0.7 size(im,1)*0.7+size(im,1)/fov*100],...
                [size(im,2)*0.9 size(im,2)*0.9],...
                'color',[1 1 1],'linewidth',5)
            %
            %             [masks,xpos,ypos,masknum] = fetchn(olf.MaskTraces & obj & 'masknum>0',...
            %                 'mask','xpos','ypos','masknum');
            %
            %             for imask = 1:length(masks)
            % %                 B = bwboundaries(masks{imask});
            % %                 plot(B{1}(:,2),B{1}(:,1),'r')
            %                  text(xpos,ypos,num2str(masknum(imask)))
            %             end
            
        end
        
        function plotTraces(obj,xl)
            
            keys = fetch(obj);
            for session = keys
                
                % compute stimuli
                figure
                hold on
                [strials,stimuli] = fetchn(olf.StimPeriods & obj,'trial','stimulus');
                [fps,ftrials] = fetch1(obj,'fps','trials');
                [sindx,desc] = fetchn(olf.Stimuli & obj,'stimulus_index','stimulus_description');
                % compute traces
                traces = fetchn(olf.Traces & session & 'trace_opt = 3','trace');
                T = cell2mat(traces');
                
                tr = double(T);
                k = hamming(round(fps/2)*1+1);
                k = k/sum(k);
                tr = convmirr(tr,k);
                
                tr2 = tr(1:1:end,:);
                fps2 = fps;
                
                % plot stimuli
                ustim = unique(stimuli);
                colors = hsv(length(ustim));
                colors = cbrewer('qual','Pastel1',length(ustim));
%                 colors(colors==0) = 0.8;
                for itrial = 1:length(strials)
                    start = find(ftrials==strials(itrial),1,'first');
                    stop = find(ftrials==strials(itrial),1,'last');
                    color = colors(strcmp(ustim,stimuli(itrial)),:);
                    area([start stop]/fps,[size(traces,1)+1 size(traces,1)+1],...
                        'facecolor',color,'edgecolor','none')
                end
                
                % plot Traces
                plot(0:1/(fps2):size(tr2,1)/(fps2) - 1/(fps2),bsxfun(@plus,tr2/2,1:size(tr2,2)),'color',[0.4 0.4 0.4])
                
                xlim([0 size(tr2,1)/fps2-1/fps2])
                ylim([0 size(tr2,2)+1])
                ylabel('Cell #')
                xlabel('Time (sec)')
                
                plot([size(tr2,1)/fps2*0.97 size(tr2,1)/fps2*0.97],[2 4],'r','linewidth',2)
                
                l = legend(desc(cellfun(@str2num,ustim)));
            end
            
            
        end
    end
end