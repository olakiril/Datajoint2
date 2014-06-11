function out = getImStats(obj,key,functions,varargin)

params.bin = 100;
params.pixels = 20;

params = getParams(params,varargin);

% get RFs
[snr gfit] = fetchn(RFFit(key),'snr','gauss_fit');

% get only good rfs and combine
gfit = mean(cell2mat(gfit(snr>mean(snr)/2)'),2);

% get the rf info
im = fetchn(RFMap(key),'on_rf');
imsize = size(im{1});

% gaussfit
rfline = gaussf(gfit);

% stim times and types of movies shown
movieNum = unique(fetchn(StatsPresents(key),'movie_num'));
if datenum(key.exp_date,'yyyy-mm-dd')>datenum('2012-02-12','yyyy-mm-dd')
    path = getLocalPath('/lab/users/philipp/stimuli/MouseMovie/MacNew/mov1l11s');
else
    path = getLocalPath('/lab/users/philipp/stimuli/MouseMovie/Mac/mov');
end

for ivar = 1:length(functions)
    eval(['val' num2str(ivar) '=cell(length(movieNum),1);']);
end

% loop through different movies shown
for iMov = 1:length(movieNum)
    display(num2str(movieNum(iMov))) 
    if strcmp(key.movie_type,'phase')
        type = 'phs';
    else
        type = 'nat';
    end
    movie = mmreader([path num2str(movieNum(iMov)) '_' type '.avi']); %#ok<TNMLP>
    nFrames = movie.NumberOfFrames;
    
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
    
    stepframes = max([1 round(movie.FrameRate*(params.bin/1000))]);
    
    % do it frame by frame
    for iFrame = 1 : stepframes: nFrames
        % get movie and select relevant pixels
        steps = iFrame:iFrame+stepframes-1;
        steps(steps>nFrames) = [];
        mov = nan(movie.Height,movie.Width,length(steps));
        for istep = 1:length(steps)
            mov(:,:,istep) = mean(read(movie,steps(istep)),3);
        end
        mov = mean(mov,3);
        mov = mov(round(ymin:ymax),round(xmin:xmax));
        mov = imresize(mov,sqrt(numel(mov)\params.pixels));
        % do the computations
        for ivar = 1:length(functions)
            eval(['val' num2str(ivar) '{iMov}(iFrame) = ' num2str(feval(functions{ivar},mov(:))) ';']);
        end
    end
end

out = cell(length(varargin),1);
for ivar = 1:length(varargin)
    eval(['out{ivar} = val' num2str(ivar) ';'])
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



