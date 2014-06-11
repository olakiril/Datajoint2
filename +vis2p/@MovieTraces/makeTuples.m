function makeTuples( obj, key )


% get trace
trace = fetch1(Traces(key),'trace');
fps    = fetch1( Movies(key), 'fps' );

% Load times of the trace
times = fetch1(VisStims(key),'frame_timestamps');

%stim times and types
stimTimes = fetchn(MoviePresents(key),'movie_times');
movieTypes = fetchn(MoviePresents(key),'movie_start_time');
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

% make sure there are no edge effects by removing 0.5 sec from the start and
% 0.5 sec from the end of the traces
traces = traces(round(fps/2):end - round(fps/2),:);

% calculate mean
key.mean_fr = mean(trace);

% loop through different bins and insert
for undersample = 0:1
    for binsize = [100 500]
        tuple = key;
        tuple.binsize = binsize;
        %bin
        d = max(1,round(tuple.binsize/1000*fps));
        k = ones(d,1)/d;
        trace = conv2(traces,k,'valid');

        tuple.undersample = undersample;
        if undersample
            trace = trace(1:d:end,:);
        end

        % get the zscore
        traceZ = zscore(trace);

        % calculate all combinations of correlations
        if length(movieTypes)>1
            cmbIn = cell(size(movieTypes,2),1);
            cmbOut = cell(size(movieTypes,2),1);
            for iMovie = 1:size(movieTypes,2)
                % get the index of the same movie trials
                trialIn = find(movieTypes(:,iMovie));
                trialOut = find(~movieTypes(:,iMovie));

                % compute permuted correlations
                cmbIn{iMovie} = combnk (trialIn,2);
                cmbOut{iMovie} = setxor(setxor(combnk([trialIn;trialOut],2), ...
                    cmbIn{iMovie},'rows'),combnk(trialOut,2),'rows');
                randOut = randperm(size(cmbOut{iMovie},1));
                cmbOut{iMovie} = cmbOut{iMovie}(randOut(1:size(cmbIn{iMovie},1)),:);
            end
            cmbIn = cell2mat(cmbIn);
            cmbOut = cell2mat(cmbOut);

            % Correlate and calculate significance
            [tuple.inCorr tuple.inCorrP] = corr(reshape(traceZ(:,cmbIn(:,1)),[],1),...
                reshape(traceZ(:,cmbIn(:,2)),[],1));
            [tuple.outCorr tuple.outCorrP] = corr(reshape(traceZ(:,cmbOut(:,1)),[],1),...
                reshape(traceZ(:,cmbOut(:,2)),[],1));
        else
            tuple.inCorr = 0;
            tuple.outCorr = 0;
            tuple.inCorrP = 0;
            tuple.outCorrP = 0;
        end
        
        % reshape trace
        trace = reshape(trace,[],1);

        % calculate CV and FF
        tuple.cv = std(trace(:))/mean(trace(:));
        tuple.ff = std(trace(:))^2/mean(trace(:));

        % calculate Lifetime sparseness
        L = length(trace);
        tuple.life_sparse = (1 - ((sum(trace)/L)^2/sum(trace.^2/L)))/(1 - (1/L));

        % calculate autocorrelation
        [auto_corr.corr auto_corr.lag auto_corr.confidence] = autocorr(trace);
        tuple.auto_corr = auto_corr;

        % calculate kurtosis
        tuple.kurtosis = kurtosis(trace);

        % insert data
        insert( obj, tuple );

    end
end



