%%
% setup params
bin = 1000;
no_dims = 3;
initial_dims = 50;
perplexity = 30;

% get data
[Data, Trials] = getData(mov3d.DecodeTime,k,bin);

% arrange data
Labels = nan(size(Trials));
Labels(1,:) = 1;
Labels(2,:) = 2;
Labels = Labels(:);
trials = Trials(:);
uTrials = unique(trials);
NewData = Data(:,:);

% prefilter data
for itrial = 1:length(uTrials);
    trialIdx = trials==uTrials(itrial);
    for icell = 1:size(Data,1)
        NewData(icell,trialIdx) = conv(NewData(icell,trialIdx),gausswin(5),'same');
    end
end

% do it
mappedX = tsne(NewData', [], no_dims, initial_dims, perplexity);

%% Plot
tbin = 3;
uTrials = unique(trials);
figure
hold on;
colors = [0 0 0.5;0.5 0 0];

for itrial = 2:100%length(uTrials);
    
    
    trialIdx = trials==uTrials(itrial);
    tIdx = all(Labels(trialIdx)==2);
     
    preIdx = trials==(uTrials(itrial)-1);
    pIdx = all(Labels(preIdx)==2);
    
    if tIdx==pIdx
        cis = 1;
    else
        cis = 0;
    end
    
    cidx = tIdx*2+cis;
 
    col = linspace(cidx+0.01,1+cidx,length(xx));
    
    xx = interpn(mappedX(trialIdx,1),tbin,'cubic');
    yy = interpn(mappedX(trialIdx,2),tbin,'cubic');
    zz = interpn(mappedX(trialIdx,3),tbin,'cubic');
    %     plot3(xx,yy,zz,'color',colors(cidx,:))
    surface([xx; xx],[yy; yy],[zz; zz],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',1);
    
end



colors = linspace(0,1,100)';
zero = ones(size(colors));
colormap([colors colors zero;colors zero zero;zero colors colors;zero zero colors]);

%% Plot
tbin = 3;
uTrials = unique(trials);
figure
hold on;
colors = [0 0 1;0.2 0.8 1;1 0 0;1 0.8 0.2];
mxx = [];
myy = [];
mzz = [];
mxx2 = [];
myy2 = [];
mzz2 = [];
for itrial = 2:length(uTrials);
    
    
    trialIdx = trials==uTrials(itrial);
    tIdx = all(Labels(trialIdx)==2);
     
    preIdx = trials==(uTrials(itrial)-1);
    pIdx = all(Labels(preIdx)==2);
    
    if tIdx==pIdx
        cis = 1;
    else
        cis = 0;
    end
    
    cidx = tIdx*2+cis;
    
    col = linspace(cidx+0.01,1+cidx,length(xx));
    
    xx = interpn(mappedX(trialIdx,1),tbin,'cubic');
    yy = interpn(mappedX(trialIdx,2),tbin,'cubic');
    zz = interpn(mappedX(trialIdx,3),tbin,'cubic');

%    p  = patchline(xx,yy,zz,'linestyle','-','edgecolor',colors(cidx+1,:),'linewidth',2,'edgealpha',0.2);
%     if tIdx
%          mxx2 = mxx2+xx;
%          myy2 = myy2+yy;
%          mzz2 = mzz2+zz;
%     else
%          mxx = mxx+xx;
%          myy = myy+yy;
%          mzz = mzz+zz;
%     end
    
          surface([xx; xx],[yy; yy],[zz; zz],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',1,'facealpha',0.5,'edgealpha',0.5);
    

%     alpha(.7)
    plot3(xx(1),yy(1),zz(1),'.','color',colors(cidx+1,:),'markersize',10)
    
end

coloridx = linspace(1,0.5,100)';

colormap([bsxfun(@times,(1-colors(1,:)),coloridx) + colors(1,:);...
    bsxfun(@times,(1-colors(2,:)),coloridx) + colors(2,:);...
    bsxfun(@times,(1-colors(3,:)),coloridx) + colors(3,:);...
    bsxfun(@times,(1-colors(4,:)),coloridx) + colors(4,:)])
 



