function makeTuples( obj, key )

tf = 0.16;

% get Optical data
disp 'loading movie...'

tpr = tpReader(Scans(key));
data = getImagingChannel(tpr,1);
data = permute(data(:,:,:),[3 1 2]);
Fs = getFramerate( tpr );

% get the vessel image
vessels = getImagingChannel(tpr,2);
vessels = mean(vessels(:,:,:),3);

imsize = size(data,2);
data = (bsxfun(@minus,double(data(:,:)),mean(double(data(:,:))))); % subtract mean for the fft

T = 1/Fs; % frame length
L = size(data,1); % Length of signal
t = (0:L-1)*T; % time series

% do it
disp 'computing...'
R = exp(2*pi*1i*t*tf)*data;
imP = squeeze(reshape((angle(R)),imsize,imsize));
imA = squeeze(reshape((abs(R)),imsize,imsize));

% save the data
disp 'inserting data...'
key.amp = imA;
key.ang = imP;
if ~isempty(vessels); key.vessels = vessels; end

insert(obj,key)

disp 'done!'

