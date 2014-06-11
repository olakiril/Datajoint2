function makeTuples( obj, key )

[decoder,th_nonlin,ntrials,nrmz] = fetch1(StatsSimParams(key),'decoder','th_nonlin','ntrials','normalize');

traces = fetch1(StatsSimTraces(key),'sim_traces');

if strcmp(th_nonlin,'rectify')
    traces =  max(0,real(traces));
elseif strcmp(th_nonlin,'absolute')
    traces =  abs(real(traces));
end

% replace nans with 0;
traces(isnan(traces)) = 0;

% normalize resposes of the filters
if nrmz;tracesp = normalize(traces,2);end

% get mean,variance,kurtosis,eu dist
key.mean = nanmean(traces,2);
key.variance = var(traces,[],2);
key.p_spars = sparseness(tracesp)';
key.l_spars = sparseness(traces')';
key.e_dist = pdist(traces')';
c = corr(traces');
key.s_corr = c(logical(tril(ones(size(c)),-1)));

% create poisson spikes
traces = poissrnd(repmat(normalize(traces),[1 1 ntrials])); %#ok<NASGU>

% classification performance
eval(['key.d_perf = ' decoder '(traces,''cells'',size(traces,1),''avg'',0);']);

% insert data
insert( obj, key );









