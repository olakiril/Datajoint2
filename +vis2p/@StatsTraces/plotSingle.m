function plotSingle(obj,cells)

% function plotSingle(obj,cells)
%
% Plots the traces of many cells

figure

for k = 1:length(cells)
    
    % get trace
    key = fetch(Traces(['masknum = ' num2str(cells(k))]).*obj);
    trace = fetch1(Traces(key),'trace');
    key = fetch(Movies.*obj);
    fps    = fetch1( Movies(key), 'fps' );
    
    % Load times of the trace
    key = fetch(VisStims.*obj);
    times = fetch1(VisStims(key),'frame_timestamps');
    
    %stim times and types
    [stimTimes movieTypes] = fetchn(StatsPresents.*obj,'movie_times','movie_num');
    movieTypes = bsxfun(@eq,movieTypes,unique(movieTypes)');
    
    % find trace segments
    traces = cell(1,length(stimTimes));
    for iTimes = 1:length(stimTimes)
        traces{iTimes} = trace(times>stimTimes{iTimes}(1) & times<stimTimes{iTimes}(end));
    end
    
    % equalize
    traces = cellfun(@(x) x(1:min(cellfun(@length,traces))),traces,'UniformOutput',0);
    
    % collapse segments
    traces = cell2mat(traces);
    
    %plot
    tt = 1/fps * ((1:size(traces,1))-.5);
    ns = sum(movieTypes,1);
    plotpos = 1:size(movieTypes,2);
    for iplot = 1:size(movieTypes,2)
        
        %Compute correlations between trials
        [i,j] = ndgrid(1:size(traces(:,movieTypes(:,iplot)),2),1:size(traces(:,movieTypes(:,iplot)),2));
        cor = corr(traces(:,movieTypes(:,iplot)));
        cor = mean(cor(i<j));
        cor = round(cor*100)/100;
        
        subplot(length(cells),size(movieTypes,2),plotpos(iplot) + size(movieTypes,2)*(k - 1))
        
        plot(tt,bsxfun(@plus,traces(:,movieTypes(:,iplot)), ...
            max(traces(:)):max(traces(:)):max(traces(:))*(ns(iplot))),'Color',[0 0 0],'linewidth',2)
        hold on
        set(gca,'Box','Off')
        set(gca,'YColor',[1 1 1])
        set(gca,'XColor',[1 1 1])
        set(gca,'TickLength',[0 0])
        set(gca,'yLim',[0 max(traces(:))*(ns(iplot)+1)])
        set(gca,'XLim',[0 ceil(tt(end))])
        set(gca,'YTick',max(traces(:)):3*max(traces(:)):max(traces(:))*(ns(iplot)))
        set(gcf,'Color',[1 1 1]);
        
        text(ceil(tt(end)),double(max(traces(:))*(ns(iplot)+1.5)),['cor:' num2str(cor)], ...
            'Color',[0.8 0.3 0.3],'HorizontalAlignment','right',...
            'FontSize',22);
        set(gca,'FontSize',22)
        
        if k == 1
            title(['Movie ' num2str(iplot)],'Color',[0 0 0],'FontSize',26);
        end
        
        if iplot == 1
            ylabel(['Cell ' num2str(k)])
            set(gca,'YColor',[0 0 0],'FontSize',26)
            plot([0 0],[0 ns(iplot)],'Color',[1 1 1],'linewidth',4)
        end
        
        pos = get(gca,'position');
        if k == 1 && iplot == 1
            set(gca,'position',[pos(1)*1.05 pos(2)*0.92 pos(3)*1.15 pos(4)*1.15])
        elseif k == 1 && iplot == 2
            set(gca,'position',[pos(1)*0.95 pos(2)*0.92 pos(3)*1.15 pos(4)*1.15])
        elseif k == 2 && iplot == 1
            set(gca,'position',[pos(1)*1.05 pos(2)*1.08 pos(3)*1.15 pos(4)*1.15])
        elseif k == 2 && iplot == 2
            set(gca,'position',[pos(1)*0.95 pos(2)*1.08 pos(3)*1.15 pos(4)*1.15])
        end
        
        if k == length(cells)
            set(gca,'XColor',[0 0 0])
            set(gca,'Xtick',1:2:9)
            xlabel('Time (sec)','FontSize',22);
            
            
            if iplot == size(movieTypes,2)
                set(gca,'YTickLabel',1:3:ns(iplot))
                set(gca,'YColor',[0 0 0])
                set(gca,'Yaxislocation','right')
                ylabel('Trials')
                set(gca,'YColor',[0 0 0],'FontSize',22)
            else
                 set(gca,'YTick',[])
            end
        else
            set(gca,'YTick',[])
        end
        
    end
end
