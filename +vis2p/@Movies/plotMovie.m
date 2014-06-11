function plotMovie(obj,tracetype)
global im
global ih
global im2
global ih2
    
if nargin<2 || isempty(tracetype)
    tracetype = 'calcium_trace';
end

%get key
key = fetch(obj);
% assert(length(key)~=1)
isAOD = strcmp(fetch1(Scans(key),'scan_prog'),'AOD');

if ~isAOD
    tp = tpReader(Scans(key));
    im = getAlignedImagingChannel(tp,1);
    im2 = getImagingChannel(tp,1);
end

if strcmp(tracetype,'norm')
    gtraces = cell2mat(fetchn(MaskTraces(key,'masknum > 0'), 'calcium_trace')');
    rtraces = cell2mat(fetchn(MaskTraces(key,'masknum > 0'),'annulus_trace')');
    traces = gtraces- rtraces;
else
    % get trace
    traces = fetchn(MaskTraces(key,'masknum > 0'),tracetype);
    traces2 = fetchn(Traces(key,'trace_opt = 17 and masknum > 0'),'trace');
    traces = cell2mat(traces');
    traces2 = cell2mat(traces2');
end

[fps, x, y]   = fetchn( Movies(key), 'fps','alignXshifts','alignYshifts');
try
    [ball, eye] = fetch1(ExtraTraces(key),'speed_trace','eye_trace');
catch
    ball = [];
    eye = [];
end
traces = calcDFoF(double(traces),fps);
scale = 1;
numpoints = 1:0.5:ceil(size(traces,2)/2)+1;
figure
 subplot(3,3,[1 2 4 5 7 8])
plot(bsxfun(@plus,traces2*scale,numpoints(1:size(traces2,2))),'k');hold on
plot(bsxfun(@plus,traces*scale,numpoints(1:size(traces,2))));
set(gca,'XTick',0:fps*100:size(traces,1),'XTickLabel',0:100:size(traces,1)/fps,...
    'YTick',2.5:2:size(traces,2)/2,'YTickLabel',4:4:size(traces,2))
plot(normalize(x{1})-2)
plot(normalize(y{1})-4,'r')
plot(ball-6,'g')
% plot(normalize(eye)-7,'y')
xlabel('time(sec)')
ylabel('cell #')
set(gcf,'Color',[1 1 1])
xlim([0 size(traces,1)]);


if ~isAOD
    % hax = axes('Units','normalized');
    pos = get(gca,'position');
    hSlider = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(traces,1),'Value',1,'Units','normalized',...
        'Position', [pos(1)-0.02 0 pos(3)+0.04 0.05],...
        'Callback', {@surfzlim,gca},'sliderstep',[0.0001 0.1]);
    addlistener(hSlider,'Value','PostSet',@(s,e) surfzlim(hSlider));


    subplot(3,3,[3])
    ih = imagesc(im(:,:,1));
    axis image
    colormap gray

    subplot(3,3,9)
    ih2 = imagesc(im2(:,:,1));
    axis image
    colormap gray
end

function surfzlim(hObj,event,ax) 
    global c

    subplot(3,3,[1 2 4 5 7 8])
    val = get(hObj,'Value');
    try
        delete(c)
    end
    c = plot([val val],get(gca,'ylim'));
    subplot(3,3,3)
    try
        delete(ih)
    end
    ih = imagesc(im(:,:,round(val)));
    axis image
     subplot(3,3,9)
    try
        delete(ih2)
    end
    ih2 = imagesc(im2(:,:,round(val)));
    axis image
end

end
   