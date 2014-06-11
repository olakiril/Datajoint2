function imTrials(obj,varargin)

params.binsize = 100; % msec
params.thr = 3; % sd's

params = getParams(params,varargin);

%get key
keys = fetch(obj);

for k = 1:length(keys)
    figure;
     set(gcf,'Position',[1 1 755 272]);
    key = keys(k);
    
    % get trace
    traces = fetchn(Traces(key),'trace');
    tracesM = cell2mat(traces');
    fps    = fetch1( Movies(key), 'fps' );
    
    % Load times of the trace
    times = fetch1(VisStims(key),'frame_timestamps');
    movietypes = [{'natural'} {'phase'}];
    
    for iMovieType = 1:2
        
        key.movie_type = movietypes{iMovieType};
        movies = fetchn(StatsPresents(key),'movie_num');
        uniMovies = unique(movies');
        
        %stim times and types
        stimTimes = fetchn(StatsPresents(key),'movie_times');
        
        % find trace segments
        traces = cell(1,length(stimTimes));
        for iTimes = 1:length(stimTimes)
            traces{iTimes} = tracesM(times>stimTimes{iTimes}(1) & times<stimTimes{iTimes}(end),:);
        end
        
        % equalize
        tracessize = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
        indx = tracessize(:,1) > mean(tracessize(:,1)) -  std(tracessize(:,1));
        traces = traces(indx);
        tracessize = tracessize(indx,1);
        traces = cellfun(@(x) x(1:min(tracessize(:,1)),:),traces,'UniformOutput',0);
        Trace = cell(1,length(uniMovies));
        for iMovie = 1:length(uniMovies)
            t = traces(movies(indx,1) == uniMovies(iMovie));
            Trace{iMovie} = cat(3,t{:});
        end
        c = cellfun(@size,Trace,'UniformOutput',0);
        if  mod(length([c{:}]),2) || sum(cell2mat(c)) == 0;continue;end
        tracessize = cell2mat(cellfun(@size,Trace,'UniformOutput',0)');
        traces = cellfun(@(x) x(:,:,1:min(tracessize(:,3))),Trace,'UniformOutput',0);
        traces = cell2mat(permute(traces,[2 1 3]));
        d = max(1,round(params.binsize/1000*fps));
        traces = convn(traces,ones(d,1)/d,'valid');
        traces = permute(traces(1:d:end,:,:),[2 1 3]);
        binsize = d*1000/fps;
        if params.thr
            sz = size(traces);
            traces = reshape(traces,size(traces,1),[]);
            sd = std(traces,[],2);
            traces = reshape(bsxfun(@gt,traces,sd),sz);
        end
        
        traces = squeeze(mean(traces));
        
        %Compute correlations between trials
        [i,j] = ndgrid(1:size(traces,2),1:size(traces,2));
        cor = corr(traces);
        cor = roundall(mean(cor(i>j)),0.01);
        
        subplot(2,1,iMovieType)
        
        imagesc(traces');
        colorbar
        hold on
        set(gca,'Box','Off')
        set(gcf,'Color',[1 1 1]);
        text(size(traces,1),0,['cor:' num2str(cor)], ...
            'Color',[0.6 0.4 0.4],'HorizontalAlignment','right','verticalAlignment','bottom',...
            'FontSize',12);
        set(gca,'FontSize',10)
        
        if iMovieType == 1
            set(gca,'XTickLabel',[])
            ylabel('Trials - Natural','FontSize',10)
        else
            set(gca,'XTickLabel',roundall((get(gca,'Xtick')*binsize)/1000,0.1))
            xlabel('Time (sec)','FontSize',12);
            ylabel('Trials - Phase','FontSize',10)
        end
       
    end
    key %#ok<*NOPRT>
    display(['movies : ' num2str(uniMovies)])
    if k ~= length(keys)
        pause
%         clf
    end
    
end

