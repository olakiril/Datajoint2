function plot3(obj,varargin)

params.trials = 1;
params.colormap = 'jet';
params.cells = 0;
params = getParams(params,varargin);

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
    movieNames = [{'Natural'} {'Phase scr'}];
    k = rmfield(key,'movie_type');
    movies = fetchn(StatsPresents(k),'movie_num');
    
    for iMovie = unique(movies');
        
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
            traces = cat(3,traces{:});
            if params.cells
               traces = traces(:,params.cells(1):params.cells(2),:); 
            end
            if ~params.trials
                traces = permute(traces,[1 3 2]);
            end
            %plot
            tt = 1/fps * ((1:size(traces,1))-.5);
            %             ns = floor(size(movies,1)/4);
            ns = size(traces,3);
            plotpos = 1:2:2*ns;
            for iplot = 1:ns
                figure(2)
                
                subplot(ns,2,plotpos(iplot) + (iMovieType - 1))
                hold on
                clims = [0 prctile(reshape(traces(:,:,iplot),[],1),95)];
                imagesc(traces(:,:,iplot)',clims)
                set(gca,'Box','Off')
                set(gca,'TickLength',[0 0])
                set(gcf,'Color',[1 1 1]);
                %                 set(gca,'Visible','off')
                set(gca,'XColor',[1 1 1])
                set(gca,'YColor',[1 1 1])
                set(gca,'YTick',[])
                set(gca,'XTick',[])
                set(gca,'XTickLabel',tt);
                set(gca,'XLim',[0 size(traces,1)]);
                set(gca,'YLim',[0 size(traces,2)])
                eval(['colormap(' params.colormap ')'])
                pos = get(gca,'position');
                set(gca,'FontSize',12)
                
                if iMovieType == 1 && iplot == 1
                    title(movieNames{iMovieType},'Color',[0.5 0.5 0.5],'FontSize',14);
                elseif iMovieType == 2 && iplot == 1
                    title(movieNames{iMovieType},'Color',[0.5 0.5 0.5],'FontSize',14);
                end
                
                if iMovie == 1
                    set(gca,'Visible','on')
                    set(gca,'XColor',[1 1 1])
                    
                    set(gca,'YColor',[0.5 0.5 0.5])
                    
                    set(gca,'position',[pos(1)*0.1 pos(2) pos(3)*1.15 pos(4)*1.2])
                    
                else
                    set(gca,'position',[pos(1)*0.9 pos(2) pos(3)*1.15 pos(4)*1.2])
                end
                
                if iplot ==ns && iMovieType == 1
                    set(gca,'Visible','on')
                    set(gca,'YColor',[1 1 1])
                    set(gca,'YTick',[])
                    set(gca,'XColor',[0.5 0.5 0.5])
                    set(gca,'YColor',[0.5 0.5 0.5])
                    set(gca,'XTick',1:50:length(tt))
                    set(gca,'XTickLabel',round(tt(1:50:length(tt))))
                    xlabel('Time (sec)','color',[0.5 0.5 0.5],'FontSize',12);
                    if ~params.trials
                        namey = 'Trials';
                    else
                        namey = 'Cells';
                    end
                    ylabel(namey,'color',[0.5 0.5 0.5],'FontSize',12);
                    
                end
            end
            
        end
        key
        pause
        clf
    end
end


