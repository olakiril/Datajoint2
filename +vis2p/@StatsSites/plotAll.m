function plotSingle(obj,varargin)

params.trials = 1;
params.cells = [];
params.thr = 0;
params.binsize = 100;
params.mean = 0;
params.pcolor = 0;

params = getParams(params,varargin);


figure

% get trace
fps    = fetch1( Movies.*obj, 'fps' );
Traces = getTraces(obj,'bin',params.binsize);

if iscell(Traces)
    lmov = length(Traces);
else
    lmov = 1;
end

for iMovie = 1:lmov
    
    traces = permute(Traces{iMovie},[2 1 3]);
    if params.cells
        traces = traces(:,params.cells,:);
    end
    
    if params.thr
        st = std(reshape(permute(traces,[1 3 2]),[],size(traces,3)));
        indx = permute(reshape(bsxfun(@lt,reshape(permute(traces,[3 2 1]),...
            size(traces,3),[]),params.thr*st'),size(traces,3),size(traces,2),size(traces,1)),[3 2 1]);
        traces(indx) = 0;
    end
    
    if ~params.trials
        traces = permute(traces,[1 3 2]);
    end
    
    if params.mean
        traces = squeeze(mean(traces,2));
    end
    %plot
    tt = params.binsize/1000* ((1:size(traces,1))-.5);
    ns = size(traces,3);
    plotpos = 1:lmov:lmov*ns - 1;
    for iplot = 1:ns
        
        subplot(ns,lmov,plotpos(iplot) + (iMovie - 1))
        hold on
        if ~params.pcolor
            clims = [0 prctile(reshape(traces(:,:,:),[],1),99)];
            imagesc(traces(:,:,iplot)',clims)
        else
            traces(traces>prctile(traces(:),99)) = prctile(traces(:),99);
            traces(traces<prctile(traces(:),1)) = prctile(traces(:),1);
            a = pcolor(double(traces(:,:,iplot)'));
            set(a,'linestyle','none')
        end
        %             set(gca,'Box','Off')
        set(gca,'TickLength',[0 0])
        set(gcf,'Color',[1 1 1]);
        %             set(gca,'XColor',[1 1 1])
        %             set(gca,'YColor',[1 1 1])
        set(gca,'YTick',[])
        set(gca,'XTick',[])
        set(gca,'XTickLabel',tt);
        set(gca,'XLim',[0 size(traces,1)]);
        set(gca,'YLim',[0 size(traces,2)])
        colormap('gray');
        pos = get(gca,'position');
        set(gca,'FontSize',12)
        
        if iMovie ~= 1 && ~params.mean
            %                 set(gca,'Visible','on')
            %                 set(gca,'XColor',[1 1 1])
            %                 set(gca,'YColor',[0.5 0.5 0.5])
            set(gca,'position',[pos(1)*0.9 pos(2) pos(3)*1.15 pos(4)*1.3])
            
        elseif ~params.mean
            set(gca,'position',[pos(1)*0.7 pos(2) pos(3)*1.15 pos(4)*1.3])
        end
        
        if iplot == 1
            title(['Movie ' num2str(iMovie)],'Fontsize',16)
        end
        set(gca,'YColor',[1 1 1])
        if iplot ==ns  && iMovie == lmov
            set(gca,'Visible','on')
            %                 set(gca,'YColor',[1 1 1])
            set(gca,'YTick',[])
            set(gca,'XColor',[0.5 0.5 0.5])
            set(gca,'YColor',[0.5 0.5 0.5])
            set(gca,'XTick',1:50:length(tt))
            set(gca,'XTickLabel',round(tt(1:50:length(tt))))
            xlabel('Time (sec)','color',[0.5 0.5 0.5],'FontSize',14);
            if ~params.trials
                namey = 'Trials';
            else
                namey = 'Cells';
            end
            ylabel(namey,'color',[0.5 0.5 0.5],'FontSize',14,'rotation',0);
            if params.mean
                namey = 'Cells';
                ylabel(namey,'color',[0.5 0.5 0.5],'FontSize',14,'rotation',90);
            end
            set(gca,'Yaxislocation','right')
            
        end
        
    end
    colormap(1- gray)
end


