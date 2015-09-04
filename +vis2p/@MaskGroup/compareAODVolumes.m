function compareAODVolumes(obj,key1,key2) %#ok<*INUSD>

% uparrow/downarrow    : move slices
% leftarrow/rightarrow : adjust AOD depth
% pageup/pagedown      : adjust ScanImage depth
% z                    : zoom in/out
% w/s/a/d              : adjust zoom window location
% escape               : quit without saving
% return               : quit and save
% spacebar             : reset offset
% o                    : set offset at current position

global data1
global data2
global olddata1
global olddata2
global pointi1 %#ok<*NUSED>
global pointi2
global offset
global markId1
global markId2
global xLoc
global yLoc
global ind
global step
global stop
global zoom
global rectH
global hardoffset
global colors
global Ids1
global Ids2

hardoffset = 0;
stop = false;
xLoc = 1;
yLoc = 1;
markId1 = [];
markId2 = [];
ind = 1;
offset = 0;
step = 1;
import vis2p.*

sigf = @(x,a,b) 1./(1+exp(-(x-b)*a));


for ikey = 1:2
    key = eval(['key' num2str(ikey)]);
    name = fetch1( Experiments*Scans(key),'file_name' );
    if ~strcmp(fetch1(Scans(key),'aim'),'stack')
        display('Detecting volume..')
        u_ind = strfind(name,'_');
        
        volumename = fetchn(Scans(rmfield(key,'scan_idx'),[...
            'file_name LIKE "' name(1:u_ind(1)) ...
            '%" and scan_idx > ' num2str(key.scan_idx -2) ...
            ' and scan_idx < ' num2str(key.scan_idx + 2) ...
            ' and aim = "Stack"' ...
            ]),'file_name');
        if isempty(volumename)
            volumename = fetchn(Scans(rmfield(key,'scan_idx'),[...
                'file_name LIKE "' name(1:u_ind(1)) ...
                '%" and scan_idx > ' num2str(key.scan_idx -3) ...
                ' and scan_idx < ' num2str(key.scan_idx + 3) ...
                ' and aim = "Stack"' ...
                ]),'file_name');
        end
        
        if isempty(volumename)
            display('No volume, skipping...')
            return
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
    [x, y, z, Ids] = fetchn(MaskCells(key),'img_x','img_y','img_z','masknum');
    points = [x, y, z];
    fn = aodReader(volumename,'Volume','ch2');
    
    % convert points to index
    pointi = nan(size(points));
    for icell = 1:size(points,1)
        [~,pointi(icell,1)] = min(abs(points(icell,1)-fn.x));
        [~,pointi(icell,2)] = min(abs(points(icell,2)- fn.y));
        [~,pointi(icell,3)] = min(abs(points(icell,3)-fn.z));
    end
    eval([ 'pointi' num2str(ikey) '=pointi;'])
    
    data = normalize(fn(:,:,:,1));
    data = repmat(reshape(data,size(data,1),size(data,2),1,size(data,3)),[1 1 3 1]);
    
    % enchance data
    data = sigf(data,3,0.1);
    data = normalize(data);
    
    eval([ 'data' num2str(ikey) '=data;'])
    
    % rectangle info
    eval(['aod_fov' num2str(ikey) '= abs(fn.x(end)-fn.x(1));'])
    [aodLoc(1), aodLoc(2)] = fetch1(Scans(key),'x','y');
    eval([ 'aodLoc' num2str(ikey) '=aodLoc;'])

    eval([ 'Ids' num2str(ikey) '=Ids;'])

end

for i = 1:size(data1,4);olddata1(:,:,:,i) = rgb2hsv(data1(:,:,:,i));end
for i = 1:size(data2,4);olddata2(:,:,:,i) = rgb2hsv(data2(:,:,:,i));end
olddata1 = data1;
olddata2 = data2;

n = round(length(pointi1)*0.8);
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

updateMasks([1 2])
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
elseif strcmp(event.Key,'z')
    zoom = ~zoom;
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
elseif C2(1,1)>0 && C2(1,1)<size(data2,1) && C2(1,2)>0 && C2(1,2)<size(data2,2)
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
global pointi1
global pointi2
global markId1
global markId2
global xLoc
global yLoc
global step
global colors
global Ids1
global Ids2

coordinates = get (gca, 'CurrentPoint');
xLoc = coordinates(1,1);
yLoc = coordinates(1,2);

pointi = eval(['pointi' num2str(event)]);
markId = eval(['markId' num2str(event)]);
data = eval(['data' num2str(event)]);
Ids = eval(['Ids' num2str(event)]);

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
    
    if event == 1; event2 = 2;else event2 = 1;end
    try emarkId2 = eval(['markId' num2str(event2) '(rindx)']);
        eval(['markId' num2str(event2) '(rindx) = [];']);
        colors(rindx) = [];end
else
    test = (sign(length(markId1)-length(markId2))+1)/2+1;
    if  test  ~= event && ~isempty(markId1) && length(markId1)~=length(markId2)
        display(['Please select a cell in channel ' num2str(test)])
        return
    end
    [~,idx] = min(pdist2([yLoc xLoc ifr],pointi));
    markId(end+1)=Ids(idx);
end


eval(['markId' num2str(event) '=markId;'])

if strcmp('alt',get(gcf,'Selectiontype'))
    updateMasks(event,emarkId,0.33,0.1)
    try    updateMasks(event2,emarkId2,0.33,0.1);end
else
    if length(markId1)==length(markId2)
        updateMasks(1,markId1(end),colors(length(markId1)))
        updateMasks(2,markId2(end),colors(length(markId1)))
        maxid = length(pointi1) + length(pointi2) - length(markId1);
        set(gcf,'name',['Found ' num2str(length(markId1)) ' common cells ' ...
           'maxID: ' num2str(maxid)])
    else
        updateMasks(event,markId(end),colors(max([length(markId1) length(markId2)])),0.3)
    end
end
updatePlot

end

function updatePlot

global ind
global offset
global data1
global data2
global pointi1
global pointi2
global markId1
global markId2
global step
global rectH
global zoom
global Ids1
global Ids2

% 
try pos = getPosition(rectH);catch; pos =[]; end

clf
for ikey = 1:2
    s = subplot(1,2,ikey);
    data = eval(['data' num2str(ikey)]);
    pointi = eval(['pointi' num2str(ikey)]);
    markId = eval(['markId' num2str(ikey)]);
    Ids = eval(['Ids' num2str(ikey)]);
    
    ifr = max(1,min(size(data,4),round(step\(ind+offset))));
    if ikey ==1;
        image(normalize(data(:,:,:,ifr)))
        celli = find(pointi(:,3)==ifr);
    else
        image(normalize(data(:,:,:,ind)))
        celli = find(pointi(:,3)==ind);
        ifr = ind;
    end
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
    set(sc,'buttondownfcn',@(s,e) eval(['cellInput(' num2str(ikey) ')']))
    title(['Slice#: ' num2str(ifr) '/' num2str(size(data,4))])
end


if ~isempty(pos)
    rectH = imrect(gca,pos);
    if zoom
        p = getPosition(rectH);
        xlim([p(1) p(1)+p(3)])
        ylim([p(2) p(2)+p(4)])
    end
end
end

function updateMasks(keys,id,Hvalue,Svalue)

global data1
global data2
global olddata1
global olddata2
global pointi1
global pointi2
global Ids1
global Ids2

if nargin<3
    Hvalue = 0.33;
    Svalue = 0.1;
    
elseif nargin<4
    Svalue = 1;
end

radius = 5 ;
for ikey = keys
    
    pointi = eval(['pointi' num2str(ikey)]);
    data = eval(['data' num2str(ikey)]);
    olddata = eval(['olddata' num2str(ikey)]);
    Ids = eval(['Ids' num2str(ikey)]);
    
    if nargin<2 || isempty(id);  idx = 1:length(pointi); 
    else     idx = find(Ids==id); end
    pointi = pointi(idx,:);
    
    for i = 1:size(data,4);data(:,:,:,i) = rgb2hsv(data(:,:,:,i));end
    maskH = squeeze(data(:,:,1,:));
    maskS = squeeze(data(:,:,2,:));
    [gy, gx, gz] = meshgrid(1:size(data,1),1:size(data,2),1:size(data,4));
    
    for i = 1:size(pointi,1)
        x = pointi(i,1);
        y = pointi(i,2);
        z = pointi(i,3);
        idx = ((gx(:) - x).^2 + (gy(:) - y).^2 + (gz(:) - z).^2) < radius^2;
        maskH(idx) = Hvalue;
        maskS(idx) = Svalue;
    end
    
    % restore grayscale information before updating the color channels
    if size(pointi,1)==1
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
    eval(['data' num2str(ikey) '=data;'])
%     if nargin<2 || isempty(id); 
%          eval(['olddata' num2str(ikey) '=data;'])
%     end
end

end