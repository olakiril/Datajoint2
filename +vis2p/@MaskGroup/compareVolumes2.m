function compareAODVolumes(obj,key1,key2) %#ok<*INUSD>

% uparrow/downarrow    : move slices
% leftarrow/rightarrow : adjust AOD depth
% pageup/pagedown      : adjust ScanImage depth
% r                    : set rectangle possition
% z                    : zoom in/out
% w/s/a/d              : adjust zoom window location
% escape               : quit without saving
% return               : quit and save
% spacebar             : reset offset
% o                    : set offset at current position

global data1
global data2
global olddata
global pointi %#ok<*NUSED>
global offset
global markId
global xLoc
global yLoc
global ind
global step
global stop
global zoom
global rectH
global hardoffset
global colors
global Ids

hardoffset = 0;
stop = false;
xLoc = 1;
yLoc = 1;
markId = [];
ind = 1;
offset = 0;
step = 1;
import vis2p.*

sigf = @(x,a,b) 1./(1+exp(-(x-b)*a));

%% AOD
name = fetch1( Experiments*Scans(key1),'file_name' );
if ~strcmp(fetch1(Scans(key1),'aim'),'stack')
    display('Detecting volume..')
    u_ind = strfind(name,'_');
    
    volumename = fetchn(Scans(rmfield(key1,'scan_idx'),[...
        'file_name LIKE "' name(1:u_ind(1)) ...
        '%" and scan_idx > ' num2str(key1.scan_idx -2) ...
        ' and scan_idx < ' num2str(key1.scan_idx + 2) ...
        ' and aim = "Stack"' ...
        ]),'file_name');
    if isempty(volumename)
        volumename = fetchn(Scans(rmfield(key1,'scan_idx'),[...
            'file_name LIKE "' name(1:u_ind(1)) ...
            '%" and scan_idx > ' num2str(key1.scan_idx -3) ...
            ' and scan_idx < ' num2str(key1.scan_idx + 3) ...
            ' and aim = "Stack"' ...
            ]),'file_name');
    end
    
    if isempty(volumename)
        display('No volume, skipping...')
        return
    else
        display(['Using: ' volumename{1}])
    end
    volumename = volumename{1};
else
    volumename = name;
end

% select an earlier nearby directory and load the volume
dirs = dir(['M:\Mouse\'  volumename(1:10) '*']);
dif = bsxfun(@minus,datenum(vertcat(dirs.name),'yyyy-mm-dd_HH-MM-SS'), datenum(volumename(1:19),'yyyy-mm-dd_HH-MM-SS'));
dirs = dirs(dif<0);
[~,idx] = max(dif(dif<0));
volumename = ['M:\Mouse\' dirs(idx).name '\AODAcq\' volumename];
[x, y, z, Ids] = fetchn(MaskCells(key1),'img_x','img_y','img_z','masknum');
points = [x, y, z];
fn = aodReader(volumename,'Volume','ch1');
aod_step_size = (fn.z(end) - fn.z(1))/length(fn.z);

% convert points to index
pointi = nan(size(points));
for icell = 1:size(points,1)
    [~,pointi(icell,1)] = min(abs(points(icell,1)-fn.x));
    [~,pointi(icell,2)] = min(abs(points(icell,2)- fn.y));
    [~,pointi(icell,3)] = min(abs(points(icell,3)-fn.z));
end

data1 = normalize(fn(:,:,:,1));
data1 = repmat(reshape(data1,size(data1,1),size(data1,2),1,size(data1,3)),[1 1 3 1]);

% rectangle info
aod_fov = abs(fn.x(end)-fn.x(1));
[aodLoc(1), aodLoc(2)] = fetch1(Scans(key1),'x','y');

%% MPScan
if nargin<3 || isempty(key2)
    key2 = [];
    key2.exp_date = key1.exp_date;
    [x, y, z] = fetch1(Scans(key1),'x','y','z');
    ks = fetch(Scans(key2,['scan_idx<' num2str(key1.scan_idx) ...
        ' and scan_prog <> "AOD" and aim = "stack"'...
        ' and x <' num2str(x+24 -14 ) ' and x >' num2str(x-10 -14) ...
        ' and y <' num2str(y+38 -28 ) ' and y >' num2str(y-10 -28) ...
        ' and z <' num2str(z+10 +80 ) ' and z >' num2str(z-90 +80)]));
    scans = [ks.scan_idx];
    if isempty(scans);disp 'No MpScan stack found, skipping...'; stop = true;return;end
    [~,idx] = min(key1.scan_idx - scans);
    key2.scan_idx = scans(idx);
end

tp = tpReader(Scans(key2));
gr = read(tp,1); gr= reshape(gr,size(gr,1),size(gr,2),1,size(gr,3));
rd = read(tp,2); rd= reshape(rd,size(rd,1),size(rd,2),1,size(rd,3));
galvo_step_size = abs(tp.hdr.acq.zStepSize);
step = aod_step_size/galvo_step_size;

data2 = normalize(single(rd));
data2(:,:,2,:) = normalize(single(gr));
data2(:,:,3,:) = zeros(size(gr));

data2 = data2(end:-1:1,end:-1:1,:,:);
% for i = 1:size(data1,4);olddata1(:,:,:,i) = rgb2hsv(data1(:,:,:,i));end
olddata = data1;


n = round(length(pointi)*0.8);
colors = linspace(0,1,n);
colors = colors(randperm(n));

% % rectangle info
% tp_pixelPitch = aod_fov2/size(data2,2);
% offsetX = (aodLoc1(1) - aodLoc2(1))/tp_pixelPitch;
% offsetY = (aodLoc1(2) - aodLoc2(2))/tp_pixelPitch;
% height = aod_fov1/tp_pixelPitch*1.1;
% pos = [size(data2,2)/2+offsetX-height/2 size(data2,2)/2-offsetY-height/2 height height];

%% Plot

h = figure('NumberTitle','off',...
    'Name','Match cells',...
    'KeyPressFcn',@dispkeyevent,...
    'WindowScrollWheelFcn', @doScroll);

updateMasks
updatePlot
% rectH = imrect(gca,pos);


end

function dispkeyevent(~, event)
global ind
global offset
global data1
global data2
global pointi
global markId1
global markId2
global xLoc
global yLoc
global step
global stop
global zoom
global rectH
global hardoffset

% step = size(data1,4)/size(data2,4);
if strcmp(event.Key,'uparrow')
    if ind<size(data2,4)
        ind = ind+1;
    end
elseif strcmp(event.Key,'downarrow')
    if ind>1
        ind = ind-1;
    end
elseif strcmp(event.Key,'leftarrow')
    offset = offset + 1;
elseif strcmp(event.Key,'rightarrow')
    offset = offset - 1;
elseif strcmp(event.Key,'pageup')
    if ind<size(data2,4)
        ind = ind+1;
    end
    offset = offset - 1;
elseif strcmp(event.Key,'pagedown')
    if ind>1
        ind = ind-1;
    end
    offset = offset + 1;
elseif strcmp(event.Key,'space')
    offset = hardoffset;
elseif strcmp(event.Key,'o')
    hardoffset = offset;
    display(['Offset set at: ' num2str(offset)])
elseif strcmp(event.Key,'return')
    stop = true;
    close all
    return
elseif strcmp(event.Key,'escape')
    disp 'Quiting ...'
    markId1 = [];
    markId2 = [];
    stop = true;
    close all
    data1 = [];
    data2 = [];
    return
elseif strcmp(event.Key,'r')
    zoom = false;
    try delete(rectH);end
    updatePlot
    subplot(1,2,2);
    display('Select rectangle')
    rectH = imrect(gca);
    display('Rectangle selected!')
    zoom = true;
elseif strcmp(event.Key,'z')
        zoom = ~zoom;
elseif strcmp(event.Key,'d')
    pos = getPosition(rectH);
    setPosition(rectH,[pos(1)+5 pos(2) pos(3) pos(4)]);
elseif strcmp(event.Key,'a')
    pos = getPosition(rectH);
    setPosition(rectH,[pos(1)-5 pos(2) pos(3) pos(4)]);
elseif strcmp(event.Key,'w')
    pos = getPosition(rectH);
    setPosition(rectH,[pos(1) pos(2)-5 pos(3) pos(4)]);
elseif strcmp(event.Key,'s')
    pos = getPosition(rectH);
    setPosition(rectH,[pos(1) pos(2)+5 pos(3) pos(4)]);
elseif strcmp(event.Key,'h')
    help compareAODVolumes
else
    display(['key: "' event.Key '" not assigned!'])
end

updatePlot

end

function doScroll(~,e)

global ind
global offset
global data1
global data2
global step
global rectH

try pos = getPosition(rectH);catch; pos =[]; end
if ~isempty(pos)
    X = pos(3);
    Y = pos(4);
else
    X = size(data2,1);
    Y = size(data2,2);
end

subplot(121)
C1 = get (gca, 'CurrentPoint');
subplot(122)
C2 = get (gca, 'CurrentPoint');
ifr = max(1,min(size(data1,4),round(step\(ind+offset))));

if C1(1,1)>0 && C1(1,1)<size(data1,1) && C1(1,2)>0 && C1(1,2)<size(data1,2)
    if  ifr<size(data1,4) && e.VerticalScrollCount>0
        offset = offset + 1;
    elseif ifr>1 &&  e.VerticalScrollCount<0
        offset = offset - 1;
    end
elseif C2(1,1) > 0 && C2(1,1) < X && C2(1,2) > 0 && C2(1,2) < Y
    if ind<size(data2,4) && e.VerticalScrollCount>0
        ind = ind+1;
        offset = offset - 1;
    elseif  ind>1 && e.VerticalScrollCount<0
        ind = ind-1;
        offset = offset + 1;
    end
else
    if ind<size(data2,4) && e.VerticalScrollCount>0
        ind = ind+1;
    elseif  ind>1 && e.VerticalScrollCount<0
        ind = ind-1;
    end
end

updatePlot

end

function cellInput ( event )
global ind
global offset
global data1
global data2
global pointi
global markId
global xLoc
global yLoc
global step
global colors
global Ids

coordinates = get (gca, 'CurrentPoint');
xLoc = coordinates(1,1);
yLoc = coordinates(1,2);

data = eval(['data' num2str(event)]);

if event ==1;
    ifr = max(1,min(size(data,4),round(step\(ind+offset))));
else
    ifr = ind;
end

if strcmp('alt',get(gcf,'Selectiontype'))
    [~,rmId] = min(pdist2([yLoc xLoc ifr],pointi));
    rindx = markId==Ids(rmId);
    emarkId = markId(rindx);
    markId(rindx) = [];
else
    [~,idx] = min(pdist2([yLoc xLoc ifr],pointi));
    markId(end+1)=Ids(idx);
end

if strcmp('alt',get(gcf,'Selectiontype'))
    updateMasks(emarkId)
else
    updateMasks(markId(end),0.1,0.6)
end
updatePlot

end

function updatePlot

global ind
global offset
global data1
global data2
global pointi
global markId
global step
global rectH
global zoom
global Ids

%
try pos = getPosition(rectH);catch; pos =[]; end

clf

s = subplot(1,2,1);
data = data1;
ifr = max(1,min(size(data,4),round(step\(ind+offset))));
image(normalize(data(:,:,:,ifr)))
celli = find(pointi(:,3)==ifr);
axis image
axis off
hold on

colors = ones(size(pointi,1),3);
color = colors(celli,:);
for icell = 1:length(celli)
    text(pointi(celli(icell),2),pointi(celli(icell),1),num2str(Ids(celli(icell))),...
        'color',color(icell,:))
end

sc = get(s,'children');
set(sc,'buttondownfcn',@(s,e) eval(['cellInput(1)']))
title(['AOD Slice#: ' num2str(ifr) '/' num2str(size(data,4))])

subplot(122)
image(normalize(data2(:,:,:,ind)))
axis image
axis off
title(['Galvo Slice#: ' num2str(ind) '/' num2str(size(data2,4))])


if ~isempty(pos)
    rectH = imrect(gca,pos);
    if zoom
        p = getPosition(rectH);
        xlim([p(1) p(1)+p(3)])
        ylim([p(2) p(2)+p(4)])
    end
end
end

function updateMasks(id,Hvalue,Svalue)

global data1
global data2
global olddata
global pointi
global Ids

if nargin<2
    Hvalue = 0.33;
    Svalue = 0.1;
    
elseif nargin<3
    Svalue = 1;
end

radius = 5 ;

data = data1;

if nargin<1 || isempty(id);  idx = 1:length(pointi);
else     idx = find(Ids==id); end
points = pointi(idx,:);

for i = 1:size(data,4);data(:,:,:,i) = rgb2hsv(data(:,:,:,i));end
maskH = squeeze(data(:,:,1,:));
maskS = squeeze(data(:,:,2,:));
[gy, gx, gz] = meshgrid(1:size(data,1),1:size(data,2),1:size(data,4));

for i = 1:size(points,1)
    x = points(i,1);
    y = points(i,2);
    z = points(i,3);
    idx = ((gx(:) - x).^2 + (gy(:) - y).^2 + (gz(:) - z).^2) < radius^2;
    maskH(idx) = Hvalue;
    maskS(idx) = Svalue;
end

% restore grayscale information before updating the color channels
if size(points,1)==1
    for i = 1:3
        odat = squeeze(olddata(:,:,i,:));
        dat = squeeze(data(:,:,i,:));
        dat(idx) = odat(idx);
        data(:,:,i,:) = dat;
    end
end

data(:,:,1,:) = maskH;
data(:,:,2,:) = maskS;
for i = 1:size(data,4);data(:,:,:,i) = hsv2rgb(data(:,:,:,i));end
data1=data;

end