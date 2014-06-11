function plot(obj )

%get key
keys = fetch(obj);

for k = 1:length(keys)

    key = keys(k);
    
    % get trace
    trace = fetch1(Traces(key),'trace');
    fps    = fetch1( Movies(key), 'fps' );

    % Load times of the trace
    times = fetch1(VisStims(key),'frame_timestamps');

    %stim times and types
    stimTimes = fetchn(StatsPresents(key),'movie_times');
    movieTypes = fetchn(StatsPresents(key),'movie_num');
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
    color = [.5 0 0; 0 .5 0; 0 0 .5];
    for iplot = 1:size(movieTypes,2)
        figure(1)
        subplot(1,size(movieTypes,2),iplot)
        plot(tt,bsxfun(@plus,traces(:,movieTypes(:,iplot)),0:max(traces(:)):max(traces(:))*(ns(iplot)-1)),'Color',[0.2 0.2 0.2])
        set(gca,'Box','Off')
        set(gca,'YColor',[1 1 1])
        set(gca,'XColor',[1 1 1])
        set(gca,'yLim',[0 max(traces(:))*ns(iplot)])

        figure(2)
        mtrace = mean(traces(:,movieTypes(:,iplot)),2);
        strace = std(traces(:,movieTypes(:,iplot)),[],2)/sqrt(ns(iplot));
        errorbar(tt,mtrace,strace,'color',color(iplot,:))
        set(gca,'Box','Off')
        set(gca,'xlim',[0 tt(end)],'ylim',[0 0.6])
        axis square
        hold on
    end
    
    figure(1)
    set(gcf,'Color',[1 1 1],'position',[100 400 600 400])

    figure(2)
    set(gcf,'Color',[1 1 1],'position',[750 400 400 400])
      
    
    pause
    clf(1),clf(2)
end
