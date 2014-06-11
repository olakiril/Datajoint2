function trialCorr(obj )

%get key
keys = fetch(obj);
figure
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
            if length(stimTimes)<10;continue;end
            for iTimes = 1:length(stimTimes)
                traces{iTimes} = tracesM(times>stimTimes{iTimes}(1) & times<stimTimes{iTimes}(end),:);
            end
            
            % equalize
            tracessize = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
            traces = cellfun(@(x) x(1:min(tracessize(:,1)),:),traces,'UniformOutput',0);
            
            %plot
            traces = squeeze(mean(cell2mat(cellfun(@(x) permute(x,[3 1 2]),traces,'uniformoutput',0)'),3))';
            if isempty(traces);continue;end
            
            %Compute correlations between trials
            m = mean(traces,2);
            blur = 400/fps; % set sampling rate to 400 fps;
            x = imresize(xcorr2(traces,m),blur);
            x = x(round(size(x,1)/2) - 100:round(size(x,1)/2) + 100,:)';% select -250 :250 msec window
            subplot(1,2*length(unique(movies')),sub2ind([2 length(uniMovies)],iMovieType,iMov))
            imagesc(x);
            hold on
            set(gca,'Ytick',[])
            tm = [1 size(x,2)/2 size(x,2)-1];
            set(gca,'Xtick',tm,'xticklabel',[-250 0 250])
            set(gca,'xminorgrid','on','xminortick','on');
            plot([size(x,2)/2 size(x,2)/2],[0 size(x,1)],'-.k','linewidth',2)
            set(gcf,'Color',[1 1 1]);
            set(gca,'FontSize',10)
            t = title(movieNames{iMovieType},'Color',[0.5 0.5 0.5],'FontSize',10);
            set(t,'Units','normalized')
            pos = get(t,'Position');
            set(t,'Position',[pos(1),pos(2)*1.01,pos(3)])
            set(gca,'ALimMode','manual')
            set(gca,'YColor',[0.5 0.5 0.5])
            set(gca,'XColor',[0.5 0.5 0.5])
            if iMovieType == 1
                pos = get(gca,'position');
                set(gca,'position',[pos(1) pos(2) pos(3) pos(4)])
                ylabel(['Movie ' num2str(iMov)],'Color',[0.5 0.5 0.8],'FontSize',12);
                set(gca,'xticklabel',[])
            else
                if iMov == 2
                    ylabel('Trials','Color',[0.5 0.5 0.5],'FontSize',10)
                end
                set(gca,'Yaxislocation','right')
                set(gca,'ytick',(0:2:size(traces,2))*blur);
                set(gca,'YTickLabel',get(gca,'Ytick')/blur);
                xlabel('Offset from Mean (msec)','FontSize',10);
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

