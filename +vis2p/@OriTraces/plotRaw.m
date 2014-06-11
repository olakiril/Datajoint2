function plotRaw( obj, key, varargin)

[params.resp_delay params.resp_period] = fetch1(OriTracesParams(key),'resp_delay','resp_period');

params.show_start = 1000;
params.show_end = 1000;

params = getParams(params,varargin);

% Load times of the trace
[time stim] = fetch1(VisStims(key),'frame_timestamps','stim_file');

%stim times and types
stimTimes = fetchn(OriPresents(key),'ori_times');
movieTypeB = fetchn(OriPresents(key),'ori_num');
movieTypesB = bsxfun(@eq,movieTypeB,unique(movieTypeB)');

keys = fetch(OriTraces(key));

for ikey = 1:length(keys)
    figure(100)
    clf
    set(gcf,'Color',[1 1 1])
    
    key = keys(ikey);
    % get trace
    trace = fetch1(Traces(key),'trace');
    fps    = fetch1( Movies(key), 'fps' );
    
    % find trace segments
    traces = cell(1,length(stimTimes));
    for iTimes = 1:length(stimTimes)
        traces{iTimes} = trace(time' > stimTimes{iTimes}(1) - params.show_start &...
            time' < stimTimes{iTimes}(end) + params.show_end);
    end
    
    % equalize
    %     indx = cellfun(@(x) length(x)>= ...
    %         floor(fps*mean(cellfun(@diff,stimTimes))/1000 - std(cellfun(@length,traces))),traces);
    indx = cellfun(@(x) length(x) >= max(cellfun(@length,traces))*9/10,traces);
    traces = traces(indx);
    movieTypes = movieTypesB(indx,:);
    movieType = movieTypeB(indx);
    traces = cellfun(@(x) x(1:min(cellfun(@length,traces))),traces,'UniformOutput',0);
    
    % collapse segments
    trace = cell2mat(traces);
    
    % equalize trials
    extra = find(sum(movieTypes,1)>min(sum(movieTypes)));
    indx = true(size(movieType,1),1);
    for i = 1:length(extra)
        indx(find(movieTypes(:,extra(i)) == 1,sum(movieTypes(:,extra(i)),1) - ...
            min(sum(movieTypes)),'last')) = false;
    end
    movieType = reshape(movieType(indx),[min(sum(movieTypes)) length(unique(movieType))]);
    trace = reshape(trace(:,indx)',min(sum(movieTypes)),length(unique(movieType)),[]);
    [~, indx] = sort(movieType(1,:));
    trace = trace(:,indx,:);
    
    orientations = stim.params.constants.orientation;
    mx = ceil(size(trace,3)/fps); % in sec
    
    for iOri = 1:size(trace,2)
        subplot(ceil(sqrt(size(trace,2))),ceil(sqrt(size(trace,2))),iOri)
        area([params.show_start*fps/1000 size(trace,3) - params.show_end*fps/1000],...
            [max(trace(:)) max(trace(:))],'FaceColor',[0.8 0.8 0.8],'LineStyle','none')
        hold on
        plot(squeeze(trace(:,iOri,:))')
        title(num2str(orientations(iOri)))
        set(gca,'Xtick',(0:1:mx)*fps)
        set(gca,'XtickLabel',0:1:mx)
        set(gca,'Xlim',[0 mx*fps],'ylim',[0 max(trace(:))])
        if iOri == size(trace,2)
            xlabel('Time from onset (sec)')
        end
    end
    figure(101)
    clf
    plot(OriTraces(key),'trials','dots')
    
    if length(keys) > 1 && ikey ~= length(keys)
        pause
    end
    
end



