function out = getStimulus(obj,key,varargin)

params.frames = 300;
params.rf = 1;

params = getParams(params,varargin);

% get RFs
[snr gfit] = fetchn(RFFit(key),'snr','gauss_fit');

% get only good rfs and combine
gfit = mean(cell2mat(gfit(snr>mean(snr)/2)'),2);

% get the rf info
im = fetchn(RFMap(key,'masknum = 1 and rf_opt_num = 1'),'on_rf');
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

out = cell(length(movieNum),1);
% loop through different movies shown
for iMov = 1:length(movieNum)
    display(num2str(movieNum(iMov)))
    if strcmp(key.movie_type,'phase')
        type = 'phs';
    else
        type = 'nat';
    end
    movie = mmreader([path num2str(movieNum(iMov)) '_' type '.avi']); %#ok<TNMLP>
    nFrames = min([movie.NumberOfFrames params.frames]);
    
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
    
    % do it frame by frame
    mov = nan(movie.Height,movie.Width,nFrames);
    for iFrame = 1 : nFrames
        mov(:,:,iFrame) = mean(read(movie,iFrame),3);
    end
    
    if params.rf
        % select pixels within the RF
        out{iMov} = mov(round(ymin:ymax),round(xmin:xmax),:);
    end
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



