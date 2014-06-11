function plotnum(obj,varargin)

params.trials = 1;
params.colormap = '(1 - gray)';
params.cells = 0;
params.textcolor = [0.1 0.1 0.1];
params.textsize = 12;
params.cellgap = 4;
params.restrict = [0 99];
params.units = 'a.u.';
params.thr = [];
params.timetick = 1;
params.bin = [];
params.movie = [];

params = getParams(params,varargin);

%get key
key = fetch(obj);

% get trace
traces = fetchn(Traces(key),'trace');
tracesM = cell2mat(traces');
fps    = fetch1( Movies(key), 'fps' );

% Load times of the trace
times = fetch1(VisStims(key),'frame_timestamps');
movietypes = [{'natural'} {'phase'}];
movieNames = [{'Natural'} {'Phase Scrambled'}];
k = rmfield(key,'movie_type');
movies = fetchn(StatsPresents(k),'movie_num');

t = reshape([traces{:}],[],1);
clims = [prctile(t,params.restrict(1)) prctile(t,params.restrict(2))];

if ~isempty(params.thr); clims(1) = params.thr;end
if isempty(params.movie);movs = unique(movies');else movs = params.movie;end

for iMovie = movs;
    figure
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
        if params.bin
            traces = trresize(traces,fps,params.bin);
            FPS = 1000/params.bin;
        else FPS = fps;
        end
        if length(params.cells) == 2
            traces = traces(:,params.cells(1):params.cells(2),:);
            celind = params.cells(1):params.cells(2);
        elseif length(params.cells) > 2
             traces = traces(:,params.cells,:);
            celind = params.cells;
        else  
            celind = 1:size(traces,2);
        end
        if ~params.trials
            traces = permute(traces,[1 3 2]);
        end
        %plot
        tt = 1/FPS * ((1:size(traces,1))-.5);
        %             ns = floor(size(movies,1)/4);
        
        ns = size(traces,3);
        cells = 1:params.cellgap:ns;
        for iplot = 1:ns
            
            subplot(ns,5,[1 2] + (iMovieType - 1)*3 + 5*(iplot - 1))
            hold on
            imagesc(traces(:,:,iplot)',clims)
            set(gca,'Box','Off')
            set(gca,'TickLength',[0 0])
            set(gcf,'Color',[1 1 1]);
            %                 set(gca,'Visible','off')
%             set(gca,'XColor',[1 1 1])
            set(gca,'YColor',[1 1 1])
            set(gca,'YTick',[])
            set(gca,'XTick',[])
            set(gca,'XTickLabel',tt);
            set(gca,'XLim',[0 size(traces,1)]);
            set(gca,'YLim',[0 size(traces,2)])
            eval(['colormap(' params.colormap ')'])
            pos = get(gca,'position');
            set(gca,'FontSize',params.textsize)
            
            if iMovieType == 1 && iplot == 1
                title(movieNames{iMovieType},'Color',params.textcolor,'FontSize',params.textsize+2);
            elseif iMovieType == 2 && iplot == 1
                title(movieNames{iMovieType},'Color',params.textcolor,'FontSize',params.textsize+2);
            end
            
            if iMovie == 1
                set(gca,'Visible','on')
                set(gca,'XColor',[1 1 1])
                set(gca,'YColor',params.textcolor)
                set(gca,'position',[pos(1)*0.1 pos(2) pos(3)*1.15 pos(4)*1.2])
            else
                set(gca,'position',[pos(1)*0.9 pos(2) pos(3)*1.15 pos(4)*1.2])
            end
            
            if sum(cells == iplot) > 0
                ylabel(num2str(celind(iplot)),'color',params.textcolor,'FontSize',params.textsize,...
                    'rotation',0,'horizontalalignment','right','verticalalignment','middle');
            end
            
            if iplot ==ns && iMovieType == 1
                set(gca,'Visible','on')
                set(gca,'YColor',[1 1 1])
                set(gca,'YTick',[])
                set(gca,'XColor',params.textcolor)
                set(gca,'YColor',params.textcolor)
                set(gca,'XTick',0:FPS*params.timetick:tt(end)*FPS)
                set(gca,'XTickLabel',round(get(gca,'xtick')/FPS))
                xlabel('Time (sec)','color',params.textcolor,'FontSize',params.textsize);
                if ~params.trials
                    namey = 'Trials';
                    namex = 'Cell #';
                else
                    namey = 'Cells';
                    namex = 'Trial #';
                end
                ylabel([namey ' {'],'color',params.textcolor,'FontSize',params.textsize-2,'interpreter','none',...
                    'rotation',0,'horizontalalignment','right','verticalalignment','middle');
%                 set(gca,'YColor',[1 1 1])
                [~,s] = suplabel(namex,'y');
                set(s,'color',params.textcolor,'fontsize',params.textsize+2,'position',[- 0.05 0.5 9])
            elseif iplot ==ns && iMovieType == 2
               hbar = colorbar;
               set(hbar,'location','eastOutside')
               set(hbar,'position',[0.9 0.11 0.014 0.06])
               set(get(hbar,'ylabel'),'string',params.units)
               barv = get(hbar,'ylim');
               set(hbar,'ytick',[barv(1) barv(2)])
               set(hbar,'yticklabel',[0 ,roundall(clims(2),0.01)])
            end   
        end
    end
end



