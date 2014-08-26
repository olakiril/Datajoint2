function plot(obj,type)

figure
import vis2p.*

[PotiIn, oriIn, PotiOut,pdm,fitVM,AllTraces] = fetchn(obj,...
    'PotiIn','oriIn','PotiOut','PdmIn','fitVM','oriTraces');

idx = PotiIn<0.05 & PotiOut>0.05;

oris = oriIn(idx);
% mtraces = maxTraces(idx);

if nargin<2; type = 'tunning';end
%%
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
% keys = fetch(obj);
% keys= keys(idx);
% pdm = pdm(idx);
% x = [0:45:315];
%
% for i = 1:length(keys)
%         subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)
%
%     [~,pd] = min(abs(pdm(i)*360 / 2 / pi - x));
%     [~,~,x1,x2]= circOri(1:8,pd);
%     traces = getTraces(CenterSurOri,'key',keys(i),'compute',1);
%     traces1 = squeeze(traces(1,2:end,:));
%     mtrace = nanmean(traces1,2);
%     etrace = nanstd(traces1,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'k')
%      hold on
%
%         traces2 = squeeze(traces(1+pd,2:end,:));
%     mtrace = nanmean(traces2,2);
%     etrace = nanstd(traces2,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'b')
%
%             traces3 = squeeze(mean(traces([1+x1 1+x2],2:end,:)));
%     mtrace = nanmean(traces3,2);
%     etrace = nanstd(traces3,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'r')
%
%     set(gca,'xtick',x)
%     xlim([0 315])
% end


%% plot all tunning curves

if strcmp(type,'tunning')
    
    keys = fetch(obj);
    keys= keys(idx);
    pdm = pdm(idx);
    x = [0:45:315];
    colors = parula(length(x));
    
    for i = 1:length(keys)
        subplot(ceil(sqrt(length(pdm))),ceil(sqrt(length(pdm))),i)
        
        
        traces = squeeze(AllTraces{i});
        
        [mx,imx] = max(mean(mean(traces(:,2:end,:),3),1));
        [~,ipd] = sort(min(abs([...
            circ_dist(x(imx)/180*pi,x/180*pi);...
            circ_dist(x(imx)/180*pi + pi,x/180*pi)...
            ])));

                
        traces1 = squeeze(traces(1,2:end,:));
        mtrace = nanmean(traces1,2);
        etrace = nanstd(traces1,[],2)/sqrt(size(oris{i},1));
        errorbar(x',mtrace,etrace,'r')
        hold on
        
        for iori = 1:length(ipd)
            traces2 = squeeze(traces(1+ipd(iori),2:end,:));
            mtrace = nanmean(traces2,2);
            etrace = nanstd(traces2,[],2)/sqrt(size(oris{i},1));
            errorbar(x',mtrace',etrace','color',colors(ipd(iori),:))
        end
        title(['pr: ' num2str(roundall(mean(reshape(traces(1+ipd(1:2),1+imx,:),[],1))/mx,0.01)) ...
            ' a: ' num2str(roundall(mean(reshape(traces(1+ipd(end-1:end),1+imx,:),[],1))/mx,0.01))...
            ' ' keys(i).exp_date ' ' num2str(keys(i).scan_idx) ' ' num2str(keys(i).masknum)])
        set(gca,'xtick',x,'box','off')
        xlim([-5 320])
        set(gca,'fontsize',6)
        axis off
    end
end

%% plot mean surround

% keys = fetch(obj);
% keys= keys(idx);
%
% x = [0:45:315];
%
%
% for i = 1:length(keys)
%         subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)
%
%     traces = squeeze(getTraces(CenterSurOri,'key',keys(i),'compute',1));
%
%     traces1 = squeeze(traces(1,2:end,:));
%     mtrace = nanmean(traces1,2);
%     etrace = nanstd(traces1,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace,etrace,'k')
%     hold on
%
%     traces2 = squeeze(mean(traces(2:end,2:end,:),1));
%     mtrace = nanmean(traces2,2);
%     etrace = nanstd(traces2,[],2)/sqrt(size(oris{i},1));
%     errorbar(x',mtrace',etrace','color',[1 0 0])
%
%     set(gca,'xtick',x)
%     xlim([-5 320])
%     axis off
% end

%% plot tunning curve

% keys = fetch(obj);
% keys= keys(idx);
% pdm = pdm(idx);
% x = [0:45:315];
% colors = jet(length(x));
% colors = [0 0 0;colors];
% mtrace = [];
% for i = 1:length(keys)
%     [~,ipd] = sort(min(abs([circ_dist(pdm(i),x/360*2*pi);circ_dist(pdm(i)+pi,x/360*2*pi)])));
%
%     traces = squeeze(getTraces(CenterSurOri,'key',keys(i),'compute',1));
%
%     % sort and normalize by max
%     traces = mean(traces(:,2:9,:),3);
%     [mx,imx] = max(traces(1,:));
%     traces = circshift(traces,[0 4-imx 0])/mx;
%
%     mtrace(i,1,:) = squeeze(traces(1,:));
%
%     for iori = 1:length(ipd)
%         mtrace(i,1+iori,:) = squeeze(traces(1+ipd(iori),:));
%     end
% end
%
% mtraces = squeeze(mean(mtrace));
% etraces = squeeze(nanstd(mtrace)/sqrt(size(mtrace,1)));
%
% for i = 1:size(mtraces,1)
%
%     errorbar(mtraces(i,:),etraces(i,:),'color',colors(i,:))
%     hold on
% end

%% plot 0 + 90 tunning curves

% keys = fetch(obj);
% keys= keys(idx);
% pdm = pdm(idx);
% x = 0:45:315;
% colors = [0 0 0;0 1 0;0 0.5 0;1 0 0;0.5 0 0];
% mtrace = [];
% for i = 1:length(keys)
%
%     [~,pd] = min(abs(pdm(i)*360 / 2 / pi - x));
%     [x0,x180,x90,x270]= circOri(1:8,pd);
%
%     traces = squeeze(getTraces(CenterSurOri,'key',keys(i),'compute',1));
%
%     % sort and normalize by max
%     traces = mean(traces(:,2:9,:),3);
%     [mx,imx] = max(traces(1,:));
%     traces = circshift(traces,[0 4-imx 0])/mx;
%
%     mtrace(i,1,:) = squeeze(traces(1,:));
%     mtrace(i,2,:) = squeeze(traces(1+x0,:));
%     mtrace(i,3,:) = squeeze(mean(traces(1+[x180],:),1));
%     mtrace(i,4,:) = squeeze(mean(traces(1+[x90],:),1));
%     mtrace(i,5,:) = squeeze(mean(traces(1+[x270],:),1));
% end
%
% mtraces = squeeze(mean(mtrace));
% etraces = squeeze(nanstd(mtrace)/sqrt(size(mtrace,1)));
%
% for i = 1:size(mtraces,1)
%
%     errorbar(mtraces(i,:),etraces(i,:),'color',colors(i,:))
%     hold on
% end


%% plot relative curve
if strcmp(type,'rel')
    x = [-135 -90 -45 0 45 90 135 180];
    keys = fetch(obj);
    keys= keys(idx);
    pdm = pdm(idx);
    mtrace = [];
    for i = 1:length(keys)

        traces = squeeze(AllTraces{i});

    %     sort and normalize by max
        traces = mean(traces(:,2:9,:),3);
        [mx,imx] = max(traces(1,:));
            mtrace(i,1,:) = squeeze( circshift(traces(1,:,:),[0 4-imx 0])/mx);
        traces = circshift(traces(2:end,:,:),[4-imx 4-imx 0])/mx;

        for iori = 1:size(traces,1)
            mtrace(i,1+iori,:) = traces(iori,:);
        end
    end
    trace = mtrace(:,2:end,4);
     mtraces = mean(trace);
      etraces = std(trace)/sqrt(size(trace,1));

        errorbar(x,mtraces,etraces)
    set(gca,'box','off')
    ylim([0.7 1])
    set(gca,'xtick',[-135 -90 -45 0 45 90 135 180])
    ylabel('Relative response')
    xlabel('Surround stimulus relative to preffered Orientation')
end

%% plot relative curve for all
if strcmp(type,'relAll')
    
    x = [-135 -90 -45 0 45 90 135 180];
    keys = fetch(obj);
    keys= keys(idx);
    mtraces = [];
    for i = 1:length(keys)
        
        traces = squeeze(AllTraces{i});
        
        %sort and normalize by max
        traces = (traces(:,2:9,:));
        [mx,imx] = max(mean(mean(traces(:,:,:),3),1));
        
        tr =  circshift(traces(2:end,:,:),[4-imx 4-imx 0])/mx;
        mtraces{i} = squeeze(tr(:,4,:));
        
    end
    
    colors = ['y' 'r' 'y' 'g' 'y' 'r' 'y' 'g'];
    xl = {'-135' '-90' '-45' '0' '45' '90' '135' '180'};
    for i = 1:length(mtraces)
        subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)
        mtrace = mean(mtraces{i},2);
        
        etrace = std(mtraces{i},[],2)/sqrt(size(mtraces{i},2));
        errorbar(x',mtrace,etrace)
        hold on
        for ic = 1:length(mtrace)
            plot(x(ic),mtrace(ic),'.','color',colors(ic))
        end
        set(gca,'xtick',[],'xticklabel',[],'box','off')
        xlim([-140 185])
        if i == (ceil(sqrt(length(mtraces)))^2 - ceil(sqrt(length(mtraces))) + 1)
            set(gca,'xtick',x,'xticklabel',x,'xticklabelrotation',90)
            
        end
    end
end
%% plus histogram
if strcmp(type,'relAllHist')
    
    keys = fetch(obj);
    keys= keys(idx);

    x = [-135 -90 -45 0 45 90 135 180];
    mtraces = [];
    Gtraces = [];
    for i = 1:length(keys)

        traces = squeeze(AllTraces{i});
        k = [];
        k.exp_date = keys(i).exp_date;
        k.scan_idx = keys(i).scan_idx;
        k.trace_opt = keys(i).trace_opt;
        k.center_sur_opt = keys(i).center_sur_opt;
%         TRACE = getTraces(vis2p.CenterSurOri,'key',k,'compute',1);
        [pd, poti]= fetchn(vis2p.CenterSurOri(k),'PdmIn','PotiIn');
        c = normalize(histc(pd/pi*180,0:45:315));
    %     c = normalize(squeeze(mean(mean(mean(TRACE(:,2:end,:,:),4),3),1)));
    % c = sparseness(squeeze(mean(mean(TRACE(1,2:end,:,:),4),1))');
        % sort and normalize by max
        traces = (traces(:,2:9,:));
        [mx,imx] = max(mean(traces(1,:,:),3));

        Gtraces(i,:) = circshift(c,[0 4-imx]);
        tr =  circshift(traces(2:end,:,:),[4-imx 4-imx 0])/mx;
        mtraces{i} = squeeze(tr(:,4,:));

    end

    colors = ['y' 'r' 'y' 'g' 'y' 'r' 'y' 'g'];
    xl = {'-135' '-90' '-45' '0' '45' '90' '135' '180'};
    for i = 1:length(mtraces)
        subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)
        mtrace = mean(mtraces{i},2);

        etrace = std(mtraces{i},[],2)/sqrt(size(mtraces{i},2));
        errorbar(x',mtrace,etrace)
        hold on
        for ic = 1:length(mtrace)
            plot(x(ic),mtrace(ic),'.','color',colors(ic))
        end
        plot(x,(Gtraces(i,:)),'k')
        set(gca,'xtick',[],'xticklabel',[],'box','off')
        xlim([-140 185])
        if i == (ceil(sqrt(length(mtraces)))^2 - ceil(sqrt(length(mtraces))) + 1)
            set(gca,'xtick',x,'xticklabel',x)

        end
    end
end
%% RF
% keys = fetch(obj);
% keys= keys(idx);
%
% for i = 1:length(keys)
%     subplot(ceil(sqrt(length(keys))),ceil(sqrt(length(keys))),i)
%
%    plotRF(obj,keys(i),'background','stimulus')
%
%
% end

%% plot relative curve for all
if strcmp(type,'relAllFit')
    
    x = [-135 -90 -45 0 45 90 135 180];
    keys = fetch(obj);
    keys= keys(idx);
    mtraces = [];
    for i = 1:length(keys)
        
        traces = squeeze(AllTraces{i});
        
        %sort and normalize by max
        traces = (traces(:,2:9,:));
        mtrace = mean(mean(traces(:,:,:),3),1);
        b = fitVonMises(double(mtrace),(0:45:315)/180*pi);
        mtrace = VonMis((0:45:315)/180*pi,b);
        [mx,imx] = max(mtrace);
        
        tr =  circshift(traces(2:end,:,:),[4-imx 4-imx 0])/mx;
        mtraces{i} = squeeze(tr(:,4,:));
        
    end
    
    colors = ['y' 'r' 'y' 'g' 'y' 'r' 'y' 'g'];
    xl = {'-135' '-90' '-45' '0' '45' '90' '135' '180'};
    for i = 1:length(mtraces)
        subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)
        mtrace = mean(mtraces{i},2);
        
        etrace = std(mtraces{i},[],2)/sqrt(size(mtraces{i},2));
        errorbar(x',mtrace,etrace)
        hold on
        for ic = 1:length(mtrace)
            plot(x(ic),mtrace(ic),'.','color',colors(ic))
        end
        set(gca,'xtick',[],'xticklabel',[],'box','off')
        xlim([-140 185])
        if i == (ceil(sqrt(length(mtraces)))^2 - ceil(sqrt(length(mtraces))) + 1)
            set(gca,'xtick',x,'xticklabel',x,'xticklabelrotation',90)
            
        end
    end
end

%% plot relative curve
if strcmp(type,'relFit')
    x = [-135 -90 -45 0 45 90 135 180];
    keys = fetch(obj);
    keys= keys(idx);
    pdm = pdm(idx);
    mtrace = [];
    for i = 1:length(keys)

        traces = squeeze(AllTraces{i});

    %     sort and normalize by max
        traces = mean(traces(:,2:9,:),3);
                mtraceF = mean(mean(traces(:,:,:),3),1);
        b = fitVonMises(double(mtraceF),(0:45:315)/180*pi);
        mtraceF = VonMis((0:45:315)/180*pi,b);
        [mx,imx] = max(mtraceF);
         
        mtrace(i,1,:) = squeeze( circshift(traces(1,:,:),[0 4-imx 0])/mx);
        traces = circshift(traces(2:end,:,:),[4-imx 4-imx 0])/mx;

        for iori = 1:size(traces,1)
            mtrace(i,1+iori,:) = traces(iori,:);
        end
    end
    trace = mtrace(:,2:end,4);
     mtraces = mean(trace);
      etraces = std(trace)/sqrt(size(trace,1));

        errorbar(x,mtraces,etraces)
    set(gca,'box','off')
    ylim([0.7 1])
    set(gca,'xtick',[-135 -90 -45 0 45 90 135 180])
    ylabel('Relative response')
    xlabel('Surround stimulus relative to preffered Orientation')
end

%% plot relative curve for all
if strcmp(type,'relAllFitNew')
    
    x = [-135 -90 -45 0 45 90 135 180];
    keys = fetch(obj);
    keys= keys(idx);
    mtraces = [];
    for i = 1:length(keys)
        
        traces = squeeze(AllTraces{i});
        
        %sort and normalize by max
        traces = (traces(:,2:9,:));
        mtrace = mean(mean(traces(:,:,:),3),1);
        b = fitVonMises(double(mtrace),(0:45:315)/180*pi);
        mtrace = VonMis((0:45:315)/180*pi,b);
        [mx,imx] = max(mtrace);
        traces = mean(traces(2:end,:,:),3);
        for iOri = 1:size(traces,1)
            b = fitVonMises(double(traces(iOri,:)),(0:45:315)/180*pi);
            traces(iOri,:) = VonMis((0:45:315)/180*pi,b);
        end
        tr =  circshift(traces,[4-imx 4-imx 0])/mx;
        mtraces{i} = squeeze(tr(:,4,:));
        
    end
    
    colors = ['y' 'r' 'y' 'g' 'y' 'r' 'y' 'g'];
    xl = {'-135' '-90' '-45' '0' '45' '90' '135' '180'};
    for i = 1:length(mtraces)
        subplot(ceil(sqrt(length(mtraces))),ceil(sqrt(length(mtraces))),i)
        mtrace = mean(mtraces{i},2);
        
        etrace = zeros(size(mtrace));
        errorbar(x',mtrace,etrace)
        hold on
        for ic = 1:length(mtrace)
            plot(x(ic),mtrace(ic),'.','color',colors(ic))
        end
        set(gca,'xtick',[],'xticklabel',[],'box','off')
        xlim([-140 185])
        if i == (ceil(sqrt(length(mtraces)))^2 - ceil(sqrt(length(mtraces))) + 1)
            set(gca,'xtick',x,'xticklabel',x,'xticklabelrotation',90)
            
        end
    end
end

%% plot relative curve
if strcmp(type,'relFitNew')
    x = [-135 -90 -45 0 45 90 135 180];
    keys = fetch(obj);
    keys= keys(idx);
    pdm = pdm(idx);
    mtrace = [];
    for i = 1:length(keys)
        
        traces = squeeze(AllTraces{i});
        
        %     sort and normalize by max
        traces = mean(traces(:,2:9,:),3);
        mtraceF = mean(mean(traces(:,:,:),3),1);
        b = fitVonMises(double(mtraceF),(0:45:315)/180*pi);
        mtraceF = VonMis((0:45:315)/180*pi,b);
        [mx,imx] = max(mtraceF);
        
        
        traces = mean(traces(2:end,:,:),3);
        for iOri = 1:size(traces,1)
            b = fitVonMises(double(traces(iOri,:)),(0:45:315)/180*pi);
            traces(iOri,:) = VonMis((0:45:315)/180*pi,b);
        end
        tr =  circshift(traces,[4-imx 4-imx 0])/mx;
        mtrace(i,:,:) = tr;
    end
    trace = mtrace(:,:,4);
    mtraces = mean(trace);
    etraces = std(trace)/sqrt(size(trace,1));
    
    errorbar(x,mtraces,etraces)
    set(gca,'box','off')
    ylim([0.7 1])
    set(gca,'xtick',[-135 -90 -45 0 45 90 135 180])
    ylabel('Relative response')
    xlabel('Surround stimulus relative to preffered Orientation')
end


%%
if strcmp(type,'hist')
    keys = fetch(obj);
    keys= keys(idx);
    mtraces = [];
    for i = 1:length(keys)

        traces = squeeze(AllTraces{i});
        traces = (traces(:,2:9,:));
        mtrace = mean(mean(traces(:,:,:),3),1);
        b = fitVonMises(double(mtrace),(0:45:315)/180*pi);
        mtraceF = VonMis((0:45:315)/180*pi,b);
        [mx,imx] = max(mtraceF);
        traces = mean(traces(2:end,:,:),3);
%         for iOri = 1:size(traces,1)
%             b = fitVonMises(double(traces(iOri,:)),(0:45:315)/180*pi);
%             traces(iOri,:) = VonMis((0:45:315)/180*pi,b);
%         end
        tr =  circshift(traces,[4-imx 4-imx 0])/max(mtraceF);
        mtraces{i} = squeeze(tr(:,4,:));
    end
    tr = cat(2,mtraces{:});

     si = (mean(tr([4],:)) - mean(tr([2 6],:)))./(mean(tr([4],:)) + mean(tr([2 6],:)));
    hist(si,40)
    hold on
    h = findobj(gca,'Type','patch');
    set(h,'linestyle','none')
    set(gca,'box','off','xlim',[-max(get(gca,'xlim')) max(get(gca,'xlim'))])
    plot([mean(si) mean(si)],get(gca,'ylim'),'r')
    plot([0 0],get(gca,'ylim'),'k')
    xlabel('Surround Effect Index (r_0 - r_9_0/r_0 + r_9_0)')
    ylabel('Cell #')
    set(gcf,'name','Surround Effect Hist')
end

function y = VonMis(x,b)
y = b(3) * exp(b(1)*(cos((x-(b(2))))-1)) + b(4)' + b(5) * exp(b(1)*(cos((pi+x-(b(2))))-1));