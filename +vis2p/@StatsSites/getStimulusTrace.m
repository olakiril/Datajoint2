function out = getStimulusTrace( obj, key )

if nargin<2
    key = fetch(obj);
    key = key(1);
end

% Load times of the traces
[times out.mouse_id out.scan_idx] = fetch1(VisStims(key).*StatsSites,'frame_timestamps','mouse_id','scan_idx');
out.fps = round(1/mean(diff(times))*1000); % calculate fps

% stim times and types of movies shown
[stimTimes out.movieNum out.movieTypes out.repeatNum]  = ...
    fetchn(StatsPresents(key),'movie_times','movie_num','movie_type','repeat_num');

% find trace segments for each stimulus 
out.stimTrace = zeros(size(times'));
for iTimes = 1:length(stimTimes)
    out.stimTrace(times > (stimTimes{iTimes}(1)) & ...
        times < (stimTimes{iTimes}(end))) = iTimes;
end



