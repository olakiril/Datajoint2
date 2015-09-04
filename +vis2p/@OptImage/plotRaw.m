function  plotRaw(obj,varargin)

% plotRaw(obj)
%
% Plots the intrinsic imaging Raw data + stimulus aquired with Intrinsic Imager
%
% MF 2015-06

params.downsample = 1/40; % downsample rate

params = getParams(params,varargin);

import vis2p.*
figure
hold on

keys = fetch(obj);

for ikey = 1:length(keys)
    
    set(gcf,'name',['OptMap ' keys(ikey).exp_date ' ' num2str(keys(ikey).scan_idx)])
    
    [data,times,trials] = getData(OptImage,keys(ikey));
    ucond = unique([trials.condIdx]);
    colors = hsv(length(ucond));
    for i = 1:length(trials)

        h = area([trials(i).start trials(i).end  ...
            trials(i).end trials(i).start],[0 0 size(data2,2)+4 size(data2,2)+4] );
      
        set(h,'FaceColor',colors(find(ucond==trials(i).condIdx),:));
        set(h,'EdgeColor','none');
    end
        mdata = mean(data);
        data = bsxfun(@rdivide,bsxfun(@minus,data,mdata),mdata);
%%
        traces =double(data2(:));
        hp = 0.1;
traces = traces + abs(min(traces(:)))+eps;
traces = traces(:,:)./convmirr(traces(:,:),hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
traces = bsxfun(@plus,traces,abs(min(traces)))+eps;
data2 = reshape(traces,size(data2));
%%

    data2 = reshape(permute(imresize(permute(data2,[2 3 1]),params.downsample),[3 1 2]),size(data2,1),[]);

    %%
    figure
hold on
     ucond = unique([trials.condIdx]);
    colors = hsv(length(ucond));
    for i = 1:length(trials)

        h = area([trials(i).start trials(i).end  ...
            trials(i).end trials(i).start],[0 0 size(data2,2)+4 size(data2,2)+4] );
      
        set(h,'FaceColor',colors(find(ucond==trials(i).condIdx),:));
        set(h,'EdgeColor','none');
    end
    plot(times,bsxfun(@plus,data2*5,(1:size(data2,2))),'color',[0.2 0.2 0.2])
%%
    if ikey ~= length(keys)
        pause
    end
end

%%

tdata = nan(length(trials)-1,14,size(data2,2));
for i = 1:length(trials)-1
%     for id = 1:size(data2,2)
        idx = find(times>trials(i).start,1,'first');
         tdata(i,:,:) = data2(idx-3:idx+10,:);
%     end
    end
    
    


% [u,~] = eigs( cov(double(data2)),1 );
% data4 = (data2 - data2*u*u');

%%
data2 = permute(imresize(permute(data,[2 3 1]),1/4),[3 1 2]);

[c, p] = princomp(data2(:,:));
traces = p(:,3:end)*c(:,3:end)';