function plot2(obj )

%get key
keys = fetch(obj);

for k = 1:length(keys)

    key = keys(k);

    % get trace
    trace = fetch1(Traces(key),'trace');
    fps    = fetch1( Movies(key), 'fps' );

    % Load times of the trace
    times = fetch1(VisStims(key),'frame_timestamps');
    movietypes = [{'natural'} {'phase'}];
    movieNames = [{'Natural'} {'Phase scr'}];
    for iMovie = 1:2

        key.movie_type = movietypes{iMovie};

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
        plotpos = 1:2:2*size(movieTypes,2);
        for iplot = 1:size(movieTypes,2)
            figure(1)

            %Compute correlations between trials
            [i,j] = ndgrid(1:size(traces(:,movieTypes(:,iplot)),2),1:size(traces(:,movieTypes(:,iplot)),2));
            cor = corr(traces(:,movieTypes(:,iplot)));
            cor = mean(cor(i<j));
            cor = round(cor*100)/100;

            subplot(size(movieTypes,2),2,plotpos(iplot) + (iMovie - 1))
            ttr = traces(:,movieTypes(:,iplot))';
            ttr = max(ttr(:)) - ttr;
            imagesc(ttr)
            colormap gray
            set(gca,'Box','Off')
            set(gca,'YColor',[1 1 1])
            set(gca,'XColor',[1 1 1])
            set(gca,'TickLength',[0 0])
%             set(gca,'yLim',[0 max(traces(:))*(ns(iplot)+1)])
%             set(gca,'XLim',[0 ceil(tt(end))])
            set(gca,'YTick',max(traces(:)):2*max(traces(:)):max(traces(:))*(ns(iplot)))
            set(gcf,'Color',[1 1 1]);
            text(ceil(tt(end)),double(max(traces(:))*(ns(iplot)+1.5)),['cor:' num2str(cor)], ...
                'Color',[0.6 0.4 0.4],'HorizontalAlignment','right',...
                'FontSize',12);
            set(gca,'FontSize',10)
            pos = get(gca,'position');
            if iMovie == 1 && iplot == 1
                title(movieNames{iMovie},'Color',[0.5 0.5 0.5],'FontSize',14);
            elseif iMovie == 2 && iplot == 1
                title(movieNames{iMovie},'Color',[0.5 0.5 0.5],'FontSize',14);
            elseif iMovie ==2 && iplot ==2
                set(gca,'XColor',[1 1 1])
            end

            if iMovie == 1
                set(gca,'YTickLabel',1:2:ns(iplot))
                set(gca,'YColor',[0.5 0.5 0.5])
                ylabel(['Movie ' num2str(iplot)],'Color',[0.5 0.5 0.5],'FontSize',10);
                set(gca,'position',[pos(1)*1.08 pos(2) pos(3)*1.1 pos(4)*1.1])
                if iplot == 1
                    set(gca,'YTick',[])
                end
            else
                set(gca,'position',[pos(1)*0.9 pos(2) pos(3)*1.1 pos(4)*1.1])
                
            end

            if iplot == size(movieTypes,2) && iMovie ==1
                set(gca,'XColor',[0.5 0.5 0.5])
                xlabel('Time (sec)','FontSize',10);
            end
        end

    end
    key
    pause
    clf(1)
end
