function makeTuples( obj, key )
tic

rsz = 0.5;

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% get the different stimulus conditions
uniStims = stim.params.constants.location;

% Get stim data
trial_conditions = horzcat(stim.params.trials.conditions);
allconditions = horzcat(stim.params.conditions.location);
startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));
endIndx = find(strcmp(vertcat(stim.eventTypes),'endStimulus'));

% get the selected trials info
conditions = nan(size(uniStims));
for icond = 1:size(uniStims,2)
    conditions(:,icond) = find(allconditions(1,:) == uniStims(1,icond) & ...
        allconditions(2,:) == uniStims(2,icond));
end
trials = struct;
for itrial = 1:length(trial_conditions)
    trials(itrial).start =  stim.events(itrial).syncedTimes(stim.events(itrial).types == startIndx)/1000;
    trials(itrial).end =  stim.events(itrial).syncedTimes(stim.events(itrial).types == endIndx)/1000;
    trials(itrial).condIdx = find(sum(trial_conditions(itrial) == conditions));
end

% get Optical data
disp 'loading movie...'
[name date] = fetch1( Scans(key),'file_name','exp_date' );
if ~isempty(strfind(name,'.h5'));fend = [];else fend = '.h5';end
path = ['M:/IntrinsicImaging/' datestr(date, 'yymmdd')];
filename = getLocalPath([path '/' name fend]);
[X,framerate] = getOpticalData(filename);

if rsz ~=1
    X = permute(imresize(permute(X,[2 3 1]),rsz),[3 1 2]);
end

times = 0:1/framerate:size(X,1)/framerate - 1/framerate;

% get the vessel image
k = [];
k.exp_date = key.exp_date;
k.mouse_id = key.mouse_id;
vesObj = Scans(k,'stim_engine = "None" and scan_prog = "Imager"');
if ~isempty(vesObj)
    scans = fetchn(vesObj,'scan_idx');
    [~,i] = min(abs(key.scan_idx - scans)); % select the closest scan
    keys = fetch(vesObj);
    name = fetch1(Scans(keys(i)),'file_name');
    vessels = squeeze(mean(getOpticalData([path '/' name '.h5'])));
else
    vessels = [];
end

% reshape
sz = size(X);
X = reshape(X, sz(1), []);

% baseline correction via interpolation
X=X./interp1q(floor(1:framerate*60:sz(1))',X(floor(1:framerate*60:sz(1)),:),[1:sz(1)]')-1;
X(isnan(X))=0;

% exclude the times 20 s before and after stimulus
ix = times > trials(1).start-20 & times < trials(end).end + 20;
times = times(ix)';
X = X(ix,:);

disp 'constructing design matrix...'
tau = 2;  % hemodynamic response time constant
nTrials = length(trials);
G = zeros(length(times), size(uniStims,2), 'single');
for iTrial = 1:nTrials
    trial = trials(iTrial);
    onset = trial.start;  % the second flip is the start of the drift
    offset = trial.end;
    cond = trial.condIdx;
    
    ix = find(times>=onset & times < offset);
    G(ix, cond) = G(ix, cond) ...
        + 1 - exp((onset-times(ix))/tau);
    
    ix = find(times>=offset & times < offset+5*tau);
    G(ix, cond) = G(ix, cond) ...
        + (1-exp((onset-offset)/tau))*exp((offset-times(ix))/tau);
    
    ix = find(times>=onset & times < offset+5*tau);
    len = length(ix);
    if ~exist('P')
        P = nan(size(uniStims,2),len+1,size(X,2),'single');
    end
    if ~isempty(ix)
        P(cond,1:len,:)=reshape(P(cond,1:len,:),len,[])+X(ix,:);
    end
end

P=P./(nTrials/size(uniStims,2));
G = bsxfun(@minus, G, mean(G));

disp 'regressing...'
[B, R2, Fp] = regressOpt(X, G, 0);
toc

% save the data
key.spot_amp = reshape(single(B)', sz(2), sz(3), []);
key.spot_r2  = reshape(single(R2), sz(2), sz(3));
key.spot_fp  = reshape(single(Fp), sz(2), sz(3));
key.spot_psth = reshape(P,size(uniStims,2),size(P,2),sz(2),sz(3));
if ~isempty(vessels); key.vessels = imresize(vessels,rsz); end

insert(obj,key)

