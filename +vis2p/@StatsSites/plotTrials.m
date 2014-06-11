function plotTrials(obj )

%get key
keys = fetch(obj);

for k = 1:length(keys)
    
    key = keys(k);
    
    % get trace
    traces = fetchn(Traces(key),'trace');
    tracesM = cell2mat(traces');
    fps    = fetch1( Movies(key), 'fps' );
    
    % Load times of the trace
    times = fetch1(VisStims(key),'frame_timestamps');
    movietypes = [{'natural'} {'phase'}];
    movieNames = [{'Nat'} {'Phs'}];
    movies = fetchn(StatsPresents(key),'movie_num');
    uniMovies = unique(movies');
    for iMov = 1:length(uniMovies);
        iMovie = uniMovies(iMov);
        for iMovieType = 1:2
            
            key.movie_type = movietypes{iMovieType};
            
            %stim times and types
            stimTimes = fetchn(StatsPresents(key,['movie_num =' num2str(iMovie)]),'movie_times');
            
            % find trace segments
            traces = cell(1,length(stimTimes));
            for iTimes = 1:length(stimTimes)
                traces{iTimes} = tracesM(times>stimTimes{iTimes}(1) & times<stimTimes{iTimes}(end),:);
            end
            
            % equalize
            tracessize = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
            traces = cellfun(@(x) x(1:min(tracessize(:,1)),:),traces,'UniformOutput',0);
            
            %plot
            tt = 1/fps * ((1:size(traces{1},1))-.5);
            traces = squeeze(mean(cell2mat(cellfun(@(x) permute(x,[3 1 2]),traces,'uniformoutput',0)'),3))';
            
            %Compute correlations between trials
            [i,j] = ndgrid(1:size(traces,2),1:size(traces,2));
            cor = corr(traces);
            cor = roundall(mean(cor(i>j)),0.01);
            
            subplot(1,2*length(unique(movies')),sub2ind([2 length(uniMovies)],iMovieType,iMov))
            wat = max(traces(:)):max(traces(:)):max(traces(:))*(size(traces,2));
            yy = bsxfun(@plus,traces,wat);
            plot(tt,yy,'Color',[0.2 0.2 0.2],'LineWidth',2)
            hold on
            set(gca,'Box','Off')
            set(gca,'YColor',[1 1 1])
            set(gca,'XColor',[1 1 1])
            set(gca,'TickLength',[0 0])
            set(gca,'yLim',[0 max(traces(:))*(size(traces,2)+1)])
            set(gca,'XLim',[0 ceil(tt(end))])
            set(gca,'YTick',max(traces(:)):2*max(traces(:)):max(traces(:))*(size(traces,2)))
            set(gcf,'Color',[1 1 1]);
            text(ceil(tt(end)),double(max(traces(:))*(size(traces,2)+1.5)),['cor:' num2str(cor)], ...
                'Color',[0.6 0.4 0.4],'HorizontalAlignment','right',...
                'FontSize',12);
            set(gca,'FontSize',10)
            t = title(movieNames{iMovieType},'Color',[0.5 0.5 0.5],'FontSize',14);
            set(t,'Units','normalized')
            pos = get(t,'Position');
            set(t,'Position',[pos(1),pos(2)*1.01,pos(3)])
            set(gca,'ALimMode','manual')
            if iMovieType == 1
                set(gca,'YTickLabel',1:2:size(traces,2))
                set(gca,'YColor',[0.5 0.5 0.5])
                set(gca,'XColor',[0.5 0.5 0.5])
                xlabel('Time (sec)','FontSize',10);
                pos = get(gca,'position');
                set(gca,'position',[pos(1) pos(2) pos(3) pos(4)])
                if iMov == 1
                    ylabel('Cells','Color',[0.5 0.5 0.5],'FontSize',10)
                end
            else
                x = xlabel(['Movie ' num2str(iMov)],'Color',[0.5 0.5 0.8],'FontSize',12);
                pos = get(x,'position');
                set(x,'position',[0,pos(2)*1.5,pos(3)])
                set(gca,'ytick',[])
                pos = get(gca,'position');
                set(gca,'position',[pos(1)-(0.08/length(unique(movies'))) pos(2) pos(3) pos(4)])
            end
        end
    end
    key %#ok<*NOPRT>
    if k ~= length(keys)
        pause
        clf
    end
end

