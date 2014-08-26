function compareChannels(obj,key1)

global data1
global pointi
global markId
global xLoc
global yLoc
global stop

stop = false;
xLoc = 1;
yLoc = 1;
markId = [];
pointi = [];

import vis2p.*

[pointi(:,1),pointi(:,2),pointi(:,3)] = fetchn(MaskCells(key1),'img_x','img_y','img_z');


[gr, rd] = fetch1(Movies(key1),'affine_green','affine_red');

data1 = normalize(double(rd));
data1(:,:,2,:) = normalize(double(gr));
data1(:,:,3,:) = zeros(size(gr));

%% Plot

figure('NumberTitle','off',...
    'Name','Find the patched cell',...
      'KeyPressFcn',@dispkeyevent,...
    'WindowButtonDownFcn',@cellInput);

image(normalize(data1(:,:,:)))
axis image
hold on

celli = find(pointi(:,3)==0);
for icell = 1:length(celli)
    text(pointi(celli(icell),1),pointi(celli(icell),2),'+',...
        'HorizontalAlignment','center','VerticalAlignment','middle')
    text(pointi(celli(icell),1),pointi(celli(icell),2),num2str(celli(icell)),...
        'VerticalAlignment','top')
end

end

function cellInput ( ~ , ~ )

global data1
global pointi
global markId
global xLoc
global yLoc

coordinates = get (gca, 'CurrentPoint');
xLoc = coordinates(1,1);
yLoc = coordinates(1,2);
   

if  sign(xLoc)>0 && sign(yLoc)>0
    if strcmp('alt',get(gcf,'Selectiontype'))
        [~,rmId] = min(pdist2([xLoc yLoc 0],pointi));
        markId(markId==rmId) = [];
    else
        [~,markId(end+1)] = min(pdist2([xLoc yLoc 0],pointi));
    end
end
    
clf
image(normalize(data1(:,:,:)))
axis image
axis off
hold on
celli = find(pointi(:,3)==0);
colors = ones(size(pointi,1),3);
colors(markId,:) = repmat([1 0 0],length(markId),1);
color = colors(celli,:);
for icell = 1:length(celli)
    text(pointi(celli(icell),1),pointi(celli(icell),2),'+',...
        'HorizontalAlignment','center','VerticalAlignment','middle','color',color(icell,:))
    text(pointi(celli(icell),1),pointi(celli(icell),2),num2str(celli(icell)),...
        'VerticalAlignment','top','color',color(icell,:))
end

end

function dispkeyevent(~, event)

global stop 
global markId

if strcmp(event.Key,'return')
    stop = true;
    close all
    return
elseif strcmp(event.Key,'escape')
    disp 'Quiting ...'
    markId = [];
    stop = true;
    close all
    return
else
    display(event.Key)
end

end

