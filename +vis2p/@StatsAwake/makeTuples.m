function makeTuples( obj, key )

movie = {'natural','phase'};
snm = {'Quiet','Active'};
[thrb, thrw, thre,key.trace_opt,key.stats_opt] = ...
    fetch1(StatsAwakeParams(key),'ball_thr','whisker_thr','eye_thr','trace_opt','stats_opt');

for imov = 1:length(movie);
    key.movie_type = movie{imov};

    for si=1:length(snm);eval(['tr' snm{si} '=[];']);end

    traces = getTraces(StatsSites(key),'compute2',1,'collapse',1);
    if isempty(traces);return;end
    b = getExtraTraces(StatsSites(key).*Scans('state = "awake"'),...
        'speed_trace','bin',fetch1(StatsSitesParams(key),'binsize'),'compute2',1,'collapse',1);
    w = getExtraTraces(StatsSites(key).*Scans('state = "awake"'),...
        'whisker_trace','bin',fetch1(StatsSitesParams(key),'binsize'),'compute2',1,'collapse',1);
    e = getExtraTraces(StatsSites(key).*Scans('state = "awake"'),...
        'eye_trace','bin',fetch1(StatsSitesParams(key),'binsize'),'function',...
        @(x) sqrt(sum(x.^2,2))-nanmedian(sqrt(sum(x.^2,2))),'compute2',1,'collapse',1);

    % filter states
    for si = 1:length(snm); % loop through 2 states, active - quiet
        if   strcmp(snm{si},'Quiet'); idx = (abs(e)<thre & abs(w)<thrw) | (abs(e)<thre & abs(b)<thrb); 
        elseif strcmp(snm{si},'Active'); idx = abs(e)<thre & abs(w)>thrw | (abs(e)<thre & abs(b)>thrb); end
        tr = cell(size(traces,1),size(traces,2));
        for iclass = 1:size(traces,2);
            for icell = 1:size(traces,1);
                tr{icell,iclass} = squeeze(traces(icell,iclass,squeeze(idx(1,iclass,:))));
            end
        end
        eval(['tr' snm{si} ' = tr;']);
    end

    % compute measures
    for si = 1:length(snm);
        key.state = snm{si};
        eval(['key.mi = nnclassRawCell(tr' snm{si} ')'';']) 
        eval(['key.mn = nanmean(nanmean(cell2mat(tr' snm{si} ''')));'])
        eval(['key.vr = nanmean(nanmean(cellfun(@nanvar,tr' snm{si} ''')));'])
        eval(['key.ps = nanmean(sparseness(cellfun(@nanmean,tr' snm{si} ''')''));'])
        eval(['key.ls = nanmean(sparseness(cellfun(@nanmean,tr' snm{si} ''')));'])
        eval(['c = corr(cellfun(@mean, tr' snm{si} ')'');'])
        eval('key.cr = nanmean(c(logical(tril(ones(size(c)),-1))));')
        eval(['key.av_trials = mean(reshape(cellfun(@length,tr' snm{si} '),[],1));'])
        eval(['key.rl = reliability(tr' snm{si} ');'])

        % insert data
        insert( obj, key );
    end
end

