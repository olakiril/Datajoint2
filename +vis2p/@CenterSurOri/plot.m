function plot(obj)

import vis2p.*

[PotiIn, oriIn, maxTraces, PotiOut,pdm] = fetchn(obj,...
    'PotiIn','oriIn','maxTraces','PotiOut','PdmIn');

%%
idx = PotiIn<0.05 & PotiOut>0.05;

oris = oriIn(idx);
mtraces = maxTraces(idx);

% %%
% colors = hsv(length(mtraces));
% x = [-135 -90 -45 0 45 90 135 180];
% for i = 1:length(mtraces)
%     subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)    
%     mtrace = mean(mtraces{i}/max(oris{i}(:)),2);
%     etrace = std(mtraces{i}/max(oris{i}(:)),[],2)/sqrt(size(mtraces{i},2));
%     errorbar(x',mtrace,etrace,'color',colors(i,:))
%     set(gca,'xtick',x)
% end

%%

% colors = hsv(length(mtraces));
% x = [0:45:315];
% for i = 1:length(mtraces)
%     subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)    
%     mtrace = mean(oris{i}/max(oris{i}(:)),1);
%     etrace = std(oris{i}/max(oris{i}(:)),[],1)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'color',colors(i,:))
%     set(gca,'xtick',x)
%     xlim([0 315])
% end

%%
keys = fetch(obj);
keys= keys(idx);
pdm = pdm(idx);
x = [0:45:315];

for i = 1:length(keys)
        subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)    
       
    [~,pd] = min(abs(pdm(i)*360 / 2 / pi - x));
    [~,~,x1,x2]= circOri(1:8,pd);
    traces = getTraces(CenterSurOri,'key',keys(i),'compute',1);
    traces1 = squeeze(traces(1,2:end,:));
    mtrace = nanmean(traces1,2);
    etrace = nanstd(traces1,[],2)/sqrt(size(oris{i},1));
    errorbar(x',mtrace,etrace,'k')
     hold on
     
        traces2 = squeeze(traces(1+pd,2:end,:));
    mtrace = nanmean(traces2,2);
    etrace = nanstd(traces2,[],2)/sqrt(size(oris{i},1));
    errorbar(x',mtrace,etrace,'b')
    
            traces3 = squeeze(mean(traces([1+x1 1+x2],2:end,:)));
    mtrace = nanmean(traces3,2);
    etrace = nanstd(traces3,[],2)/sqrt(size(oris{i},1));
    errorbar(x',mtrace,etrace,'r')
    
    set(gca,'xtick',x)
    xlim([0 315])
end


%%
% keys = fetch(obj);
% keys= keys(idx);
% pdm = pdm(idx);
% x = [0:45:315];
% 
% for i = 1:length(keys)
%         subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)    
%        
%     [~,ipd] = sort(abs(pdm(i)*360 / 2 / pi - x));
%     
%     
%     traces = squeeze(getTraces(CenterSurOri,'key',keys(i),'compute',1));
%     traces1 = squeeze(traces(1,2:end,:));
%     mtrace = nanmean(traces1,2);
%     etrace = nanstd(traces1,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'k')
%     hold on
%      
%     traces2 = squeeze(traces(1+ipd,2:end,:));
%     mtrace = nanmean(traces2,3);
%     etrace = nanstd(traces2,[],3)/sqrt(size(oris{i},1));
%     errorbar(mtrace',etrace')
%     
%     traces3 = squeeze(mean(traces([1+x1 1+x2],2:end,:)));
%     mtrace = nanmean(traces3,2);
%     etrace = nanstd(traces3,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'r')
%     
%     set(gca,'xtick',x)
%     xlim([0 315])
% end