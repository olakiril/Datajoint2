function plotEye(obj,clim)

eye_dir = findEye(ExtraTraces,fetch(obj));
load([eye_dir 'eyeData.mat'])
load([eye_dir 'whiskerData.mat'])
xyloObj = VideoReader([eye_dir 'eyemovie.avi']);

fn = [eye_dir 'activity.h5'];
activity =loadHWS(fn,'activity','activity');

%% filter eye movements
for ia = 1:2
    x = eyeData.centers(:,ia);
    x(x ==0|isnan(x)) = nanmean(x);
    mx = mean(x);
    x = x - mx;
    
    % high amplitude noise
    iters = 200;
    thr =20; % max shift allowed between frames in pixels
    for i = 1:iters
        d = [0; diff(x)];
        indb = ((d>thr&x>thr));
        inda = d<-thr&x<-thr ;
        x(inda) =x([inda(1);inda(1:end-1)]);
        x(indb) = x([indb(2:end);indb(end)]);
    end
    
    % low amplitude noise
    thr = 3; % max shift allowed between frames in pixels
    d = [0; diff(x)];
    ind = d>thr&x>thr | d<-thr&x<-thr;
    x(ind) = interp1(find(~ind),x(~ind),find(ind),'spline');
    
    % filter
    windowSize = 20;
    xx = conv(x,ones(1,windowSize)/windowSize,'same');
    xx(1:ceil(windowSize/2)) = mean(xx(1:ceil(windowSize/2)));
    xx(end-ceil(windowSize/2):end) = mean(xx(end-ceil(windowSize/2):end));
    ec(:,ia) = xx+mx;
end

%% filter pupil dilation
x =  eyeData.radii;
x(x ==0|isnan(x)) = nanmean(x);
mx = nanmean(x);
x = x - mx;

% high amplitude noise
iters = 200;
thr =20; % max shift allowed between frames in pixels
for i = 1:iters
    d = [0; diff(x)];
    indb = ((d>thr&x>thr));
    inda = d<-thr&x<-thr ;
    x(inda) =x([inda(1);inda(1:end-1)]);
    x(indb) = x([indb(2:end);indb(end)]);
end

% low amplitude noise
thr = 3; % max shift allowed between frames in pixels
d = [0; diff(x)];
ind = d>thr&x>thr | d<-thr&x<-thr;
x(ind) = interp1(find(~ind),x(~ind),find(ind),'spline');

% filter
windowSize = 20;
xx = conv(x,ones(1,windowSize)/windowSize,'same');
xx(1:ceil(windowSize/2)) = mean(xx(1:ceil(windowSize/2)));
xx(end-ceil(windowSize/2):end) = mean(xx(end-ceil(windowSize/2):end));
er = xx+mx;

%% plot
wt = sqrt((whiskerData.shifts.^2));
% wt(1) = wt(2);
% wt = normalize(wt)*10;
wt2 = trresize(wt,60,500);
wt = normalize(interp1(wt2,0.51:length(wt2)/length(wt):length(wt2)+0.499));
ep = sqrt(sum(ec.^2,2));
ep = ep - ep(1);
la = length(activity);
f = figure(1);
set(f,'KeyPressFcn',@pausef)
global iframe
global PS
global step
global go
clf
PS = false;
iframe = 0;
step = 10;
go = true;

while iframe<xyloObj.NumberOfFrames && go
    if ~PS
        iframe = iframe+step;
    end
    im = read(xyloObj,iframe);
    clf
    subplot(5,1,[1:3])
    if nargin<2
    imagesc(im)
    else
        imagesc(im(:,:,1),clim)
        colormap gray
    end
    title([num2str(round(15*step/60)) 'X'])
    hold on
    [x, y] = circlepoints(er(iframe));
    plot(x+ec(iframe,1)+eyeData.xind, y+ec(iframe,2)+eyeData.yind, '-w','linewidth',2);
    axis image
    
    subplot(5,1,[4 5]);
%     plot(1:10:la,ep(1:10:end));
    hold on
    plot(1:10:la,activity(1:10:end)-5,'r')
        plot(1:10:la,wt(1:10:end),'g')
%         plot([0 la],[150 150],'b')
%         plot(1:100:la,er(1:100:end),'m')
    set(gca,'buttonDownFcn',@mousepos)
%     set(gca,'xlim',[0 la])
    plot([iframe iframe],get(gca,'ylim'),'-.k')
%     legend({'Eye Position','Ball','Whiskers','Pupil'})
    drawnow
end

function mousepos(hObject,~)
global iframe
pos=get(hObject,'CurrentPoint');
iframe = round(pos(1));

function pausef(~,event)

global PS
global step
global go

if strcmp(event.Key,'space')
    PS = ~PS;
elseif strcmp(event.Key,'equal')
    step = step*2;
elseif strcmp(event.Key,'hyphen')
    if step>=2
        step = round(step/2);
    end
elseif strcmp(event.Key,'escape')
    PS = false;
    go = false;
end
