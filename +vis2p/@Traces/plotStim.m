function plotStim(obj,opts,cind,lp,fig)

if nargin<2;  opts = [3]; end
if nargin<4 || isempty(lp); lp = 10; end
% if nargin<5;  figure; end

% get traces and frames
sp = fetch1(Scans.*obj,'scan_prog');
if strcmp(sp,'AOD')
    [traces(:,:,1),traces(:,:,2)] = ...
        fetchn(MaskTraces('masknum>0').*obj,...
        'calcium_trace','red_trace');
else
    [traces(:,:,1),traces(:,:,2),traces(:,:,3)] = ...
        fetchn(MaskTraces('masknum>0').*obj,...
        'calcium_trace','red_trace','annulus_trace');
end
if nargin<3 || isempty(cind)
    cind = 1:size(traces,1);
elseif length(cind)==1; cind = 1:min([size(traces,1),cind]);
end
traces = traces(cind,:,:);
traceTypes = size(traces,3);
traces = ([traces{:}]);  % traces are in columns
traces = reshape(traces,size(traces,1),[],traceTypes);

if strcmp(sp,'AOD');traces = (traces)+abs(min(traces(:)));end

fps = fetch1(Movies.*obj,'fps');
key = fetch(Scans.*obj);
key.trace_opt = opts(1);
traceOpt = fetch(TracesOpt(key),'*');
options = fetch1(TracesOpt(key),'trace_computation');
options = strread(options, '%s','delimiter',',');

% low pass filter
traceOpt.lowPass = lp;
options{2} = 'traceFilter';

% filter traces
tracesP = traces;
for iopt = 1:length(options)
    [tracesP]= eval([options{iopt} '(tracesP,fps,traceOpt)']);
end

TRACES = tracesP;

% Load times of the trace
times = fetch1(VisStims(key).*StatsParams,'frame_timestamps');

% equalize traces 2 times if length difference is 1
if abs(length(TRACES) - length(times)) == 1
    ml = min([size(TRACES,1)  length(times)]);
    TRACES = TRACES(1:ml,:);
    times = times(1:ml);
end

% stim times and types
[stimTimes, movies, types] = fetchn(StatsPresents(key),'movie_times','movie_num','movie_type');
uniMovies = unique(movies');

% plot
colorsN = [0.6 0.6 0.8;0.7 0.7 1];
colorsP = [0.8 0.6 0.6;1 0.7 0.7];

hold on
for i = 1:length(types)
    istart = find(times>stimTimes{i}(1),1,'first');
    iend = find(times<stimTimes{i}(end),1,'last');
    h = area([istart iend  ...
        iend istart]/fps,[0 0 size(TRACES,2)+4 size(TRACES,2)+4] );
    if strcmp(types{i},'natural')
        colors = colorsN(movies(i)==uniMovies,:);
    else
        colors = colorsP(movies(i)==uniMovies,:);
    end
    set(h,'FaceColor',colors);
    set(h,'EdgeColor','none');
end

% figure
plot(0:1/fps:length(TRACES)/fps - 1/fps,bsxfun(@plus,TRACES*2,(1:size(TRACES,2))),'color',[0.2 0.2 0.2])

