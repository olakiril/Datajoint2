function makeTuples( obj, key )

key = rmfield(key,'stim_idx');

[imrz, wrf, matrix] = fetch1(StatsImagesParams(key),'imresize','wrf','matrix');

if strcmp(wrf,'yes')
    % get traces of all neurons of a site and fps
    if ~isempty(RFFit(key))
        [snr, gfit] = fetchn(RFFit(key),'snr','gauss_fit');
    else
        display('nor RF maps!')
        return
    end
    
    % get only good rfs and combine
    gfit = mean(cell2mat(gfit(snr>=mean(snr)/2)'),2);
    
    % get the rf info
    im = fetchn(RFMap(key),'on_rf');
    imsize = size(im{1});
    rfline = gaussf(gfit);
end

% stim times and types of movies shown
movieNum = unique(fetchn(StatsPresents(key),'movie_num'));
sf = fetch1(VisStims(key,'exp_type = "MoviesExperiment" or exp_type = "MouseStatExperiment"'),'stim_file');
if strcmp(sf.params.constants.movieName,'mov1l11s')
    path = getLocalPath('Q:/MouseMovie/MacNew/mov1l11s');
else
    path = getLocalPath('Q:/MouseMovie/Mac/mov');
end

% loop through different movies shown
for iMov = 1:length(movieNum)
    display(num2str(movieNum(iMov)))
    subkey = key;
    subkey = rmfield(subkey,'trace_opt');
    subkey.movie_num = movieNum(iMov);
    if strcmp(key.movie_type,'phase')
        type = 'phs';
    else
        type = 'nat';
    end
    filename = [path num2str(movieNum(iMov)) '_' type '.avi'];
    movie = VideoReader(filename); %#ok<TNMLP>
    nFrames = movie.NumberOfFrames;
    
    if strcmp(wrf,'yes')
        % get the rf possitions
        bp = rfline;
        bp(1,:) = (bp(1,:)*movie.Width)/imsize(2)+0.5;
        bp(2,:) = (bp(2,:)* movie.Height)/imsize(1)+0.5;
        maxx = max(max(bp,[],2) -  min(bp,[],2));
        xycenter = mean(bp,2);
        xmin = max(0,xycenter(1) - maxx/2);
        xmax = min(movie.Width,xycenter(1) + maxx/2);
        ymin = max(0,xycenter(2) - maxx/2);
        ymax = min(movie.Height,xycenter(2) + maxx/2);
    end
    
    % do it frame by frame
    sampleSize = round(100000*imrz);
    maxdist = round(100*imrz);
    euDist = @(x,y) sqrt(diff(x).^2 + diff(y).^2 );
    Mov = nan(movie.H,movie.W,nFrames);
    for i = 1:nFrames
        Mov(:,:,i) = mean(read(movie, i),3);
    end
    
    % get movie and select relevant pixels
    if strcmp(wrf,'yes')
        Mov = Mov(round(ymin:ymax),round(xmin:xmax),:);
    end
    if imrz~=1
        Mov = imresize(Mov,imrz,'nearest');
    end
    %     dat = permute(Mov,[3 1 2]);
    
    if strcmp(matrix,'no')
        ik = nan(nFrames,1); im = ik; is = ik;its = ik;
        ic = nan(nFrames,maxdist); id = ic;
%         dMov = diff(Mov,[],3);
%         dMov(:,:,end+1) = nan(size(Mov,1),size(Mov,2));
%         
  
        for iFrame = 1 : nFrames
            
            mov = Mov(:,:,iFrame);
%             dim = dMov(:,:,iFrame);
            % first order statistics
            im(iFrame) = mean(mov(:));
            is(iFrame) = std(mov(:));
%             mv = mov(:);
            %         mv(mv == max(mv) | mv == min(mv)) = []; % make sure hard edges don't interfere
%             ik(iFrame) = kurtosis(mv(:));
%             idd(iFrame) = mean(dim(:));
            %         % second order statistics
            %         x = randi(size(mov,1),[2 sampleSize]);
            %         y = randi(size(mov,2),[2 sampleSize]);
            %         dist = round(euDist(x,y));
            %         udist = unique(dist);
            %         for idist = 1:maxdist
            %             ic(iFrame,idist) = corr(...
            %                 mov(sub2ind(size(mov),x(1,dist == udist(idist)),y(1,dist == udist(idist))))',...
            %                 mov(sub2ind(size(mov),x(2,dist == udist(idist)),y(2,dist == udist(idist))))');
            %             id(iFrame,idist) = udist(idist);
            %         end
            %
            %         % third order statistics
            %         angles = [0 45 90];
            %         data = nan(length(angles),1);
            %         for iangle = 1:length(angles)
            %             ima = imrotate(mov,angles(iangle),'crop');
            %             bi = bispecd(reshape(ima,size(ima,1),[]));
            %             co = fftshift(real(ifftn(fftshift(bi))));
            %             imcut = round(length(co)/4);
            %             data(iangle) = mean(mean(abs(co(imcut:end - imcut,imcut:end - imcut))));
            %         end
            %         its(iFrame) = mean(data);
            %         indx = true(nFrames,1);
            %         indx(iFrame) = false;
            %         subkey.im_pdist(iFrame) = mean(pdist2(dat(iFrame,:),dat(indx,:)));
        end
        
%         if strcmp(wrf,'yes')
%             [amp, ang,pwr,mot] = ...
%                 extractMotion(filename,round([xmin ymin xmax-xmin ymax-ymin]));
%         else
%             [amp, ang,pwr,mot] = extractMotion(filename);
%         end
    else
        im = Mov(:,:,1);
        sz = eval(matrix);
        xstep = round(size(im,1)/sz(1));
        ystep = round(size(im,2)/sz(2));
        xind = 1:xstep:size(im,1);
        yind  = 1:ystep:size(im,2);
        imi = zeros(size(im));
        ind = 0;
        for iy = 1:length(yind);
            for ix = 1:length(xind);
                ind = ind+1;
                imi(xind(ix):xind(ix)+xstep,yind(iy):yind(iy)+ystep) = ind;
            end;
        end
        % correct extra point
        imi = imi(1:size(im,1),1:size(im,2));
        
        % initialize
        imis = unique(imi);
        ik = nan(length(imis),nFrames); im = ik; idd = ik;amp = ik;ang = ik;pwr = ik;
        ic = []; id = []; its = [];
        dMov = diff(Mov,[],3);
        dMov(:,:,end+1) = nan(size(Mov,1),size(Mov,2));
        
        % loop through frames
        for iFrame = 1 : nFrames
            mov = Mov(:,:,iFrame);
            dim  = dMov(:,:,iFrame);
            % loop through positions in the image
            for ipos = 1:numel(imis)
                im(ipos,iFrame) = mean(mov(imi==ipos));
                is(ipos,iFrame) = std(mov(imi==ipos));
                ik(ipos,iFrame) = kurtosis(mov(imi==ipos));
                idd(ipos,iFrame) = mean(dim(imi==ipos));
            end
        end
        if strcmp(wrf,'yes')
            [~,~,~,~,ofmov] = extractMotion(filename,round([xmin ymin xmax-xmin ymax-ymin]));
        else
            [~,~,~,~,ofmov] = extractMotion(filename);
        end

        
        for iFrame = 2 : nFrames
            for ipos = 1:numel(imis)
                of = ofmov(:,:,iFrame);
                of = of(imi==ipos);
                Xvec = nansum(real(of));       % Vectors X
                Yvec = nansum(imag(of));       % Vectors Y              
                amp(ipos,iFrame-1) = sqrt(Xvec.^2 + Yvec.^2)/numel(of);
                ang(ipos,iFrame-1) = atand(abs(Xvec/Yvec));
                pwr(ipos,iFrame-1) = nanmean(abs(of(:)));
                mot.up(iframe) =  nanmean(Xvec(Xvec<0));
                mot.down(iframe) =  nanmean(Xvec(Xvec>0));
                mot.left(iframe) =  nanmean(Yvec(Yvec<0));
                mot.right(iframe) =  nanmean(Yvec(Yvec>0));
            end
        end
    end
    
    % fill in
%     subkey.im_kurtosis = ik;
    subkey.im_mean = im;
    subkey.im_std = is;
%     subkey.im_pwz_corr = ic;
%     subkey.im_pwz_dist = id;
%     subkey.im_bispectrum = its;
%     subkey.im_motion.amp = amp;
%     subkey.im_motion.ang = ang;
%     subkey.im_motion.pwr = pwr;
%     subkey.im_motion.dir = mot;
%     subkey.im_diff = idd;          
    
    % insert data
    insert( obj, subkey );
end


function rfline = gaussf(gfit)

gm = gfit(1:2); C = diag(gfit(3:4)); cc=gfit(5)*sqrt(prod(gfit(3:4))); C(1,2)=cc; C(2,1)=cc;
npts=50;
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C);
d(d<0) = 0;
d = 2 * sqrt(d); % convert variance to sdwidth*sd
rfline = (v*d*ap) + repmat(gm, 1, size(ap,2));


