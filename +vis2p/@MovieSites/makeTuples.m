function makeTuples( obj, key )

for sig = [0.05 1]
    for qual = [4 3 2 1 0];
        % get traces of all neurons of a site and fps
        tracesR = fetchn((Traces(key).*RFStats(['onpoff_p <' num2str(sig)]))...
            .*MaskTracesQuality(['ca_event_snr>' num2str(qual)]),'trace');
        if isempty(tracesR)
            continue
        end
        tracesR = [tracesR{:}];
        fps    = fetch1( Movies(key), 'fps' );

        % Load times of the traces
        times = fetch1(VisStims(key),'frame_timestamps');

        % stim times and types of movies shown
        stimTimes = fetchn(MoviePresents(key),'movie_times');
        movieTypes = fetchn(MoviePresents(key),'movie_start_time');
        movieTypes = bsxfun(@eq,movieTypes,unique(movieTypes)');

        % find trace segments for each stimulus and remove 0.5 from each
        % side  
        traces = cell(1,length(stimTimes));
        for iTimes = 1:length(stimTimes)
            traces{iTimes} = tracesR(times > (stimTimes{iTimes}(1) + 500) & ...
                times < (stimTimes{iTimes}(end) - 500 ),:);
        end

        % remove incomplete trials
        L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
        if isempty(L)
            keyboard
        end
        indx = L(:,1)>=prctile(L(:,1),10);
        traces = traces(indx);
        L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
        traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);

        % mean across same stimulus of different repetitions and collapse segments
        traces = cat(3,traces{:});

        % insert fot each bin
        for binsize = [100 500]
            tuple =key;
            tuple.binsize = binsize;
            tuple.pThr = sig;
            tuple.snrThr = qual;

            % rebin to approximate binsize (downsampling by an integer factor)
            d = max(1,round(tuple.binsize/1000*fps));
            tuple.actual_binsize = d*1000/fps;
            k = ones(d,1)/d;
            tr = convn(traces,k,'valid');
            tr = tr(1:d:end,:,:);

            % compute cross-correlogram (before downsampling)
            c = cellfun(@corr,num2cell(tr,[1 2]),'UniformOutput',0);
            [i,j] = ndgrid(1:size(c{1},1),1:size(c{1},1));
            c = cell2mat(c);
            g = cat(1,arrayfun(@(x,y) x(y),c,repmat(i<j,[1 1 size(c,3)]),'UniformOutput',false));
            tuple.cor = squeeze(nanmean(nanmean(cellfun(@mean,g),1),2));

            %reshape
            tr = reshape(permute(tr,[1 3 2]),[],size(tr,2));
            
            % kurtosis
            tuple.kurtosis = kurtosis(tr')';

            % compute population sparseness
            pops = @(x,y) (1 - ((sum(x)/y).^2./sum(x.^2/y)))/(1 - (1/y));
            tuple.pop_sparse = pops(tr',size(tr,2))';

            % find stim length
            tuple.stim_length = size(tr,1)/size(movieTypes,2);

            % calculate avtivity sparseness
            tuple.act_sparse = mean(bsxfun(@lt,tr,mean(tr,2) + std(tr,1,2)) & bsxfun(@gt,tr,mean(tr,2) - std(tr,1,2)),2);

            % compute synchrony
            v = convmirr( (tr - convmirr(tr.^2,k)).^2, k );
            v = mean(v,2);
            m = mean(tr,2);
            mv = convmirr( (m - convmirr(m.^2,k)).^2, k );
            tuple.synchrony = mv./v;

            % insert data
            insert( obj, tuple );
        end
    end
end