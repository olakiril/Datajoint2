function makeTuples( obj, key)

tuple = key;

import vis2p.*

state = strcmp(fetch1(Scans(key),'stim_engine'),'State');
if state
    if ~isempty(Scans(key,'exp_date > "2012-09-01"'))
        stimPath = '/stor01/stimulation/Mouse';
    else
        stimPath = '/stor01/stimulation/Mouse2P1';
    end
else
    stimPath = '/stor01/stimulation/TwoPhoton';
end
scpr = fetch1(Scans(key),'scan_prog');
if strcmp(scpr,'AOD') || strcmp(scpr,'Unirec')
    [path,name] = fetch1( Experiments*Scans(key), 'directory','file_name' );
    if ~isempty(strfind(name,'.h5'));fend = [];else fend = '.h5';end
    filename = getLocalPath([path '/' name fend]);
    isSyncedTpr = 0;
    if ~isempty(Scans(key,'exp_date > "2012-05-01"'))
        dirs = dir(['M:\Mouse\' name(1:10) '*']);
        names = vertcat(dirs.name);
        % find the session that started immediatelly before the recording
        timediff = str2num(names(:,[12 13 15 16 18 19]))- str2num( name([12 13 15 16 18 19]));
        timeind = find(timediff<0);
        [foo,i] = max(timediff(timeind));
        itime = timeind(i);
        filename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' name];
        aodReader = 'new';
    end
elseif strcmp(scpr,'MPScan')
    tpr = tpReader( Scans(key) );
    filename = getLocalPath(getFilename(tpr));
    isSyncedTpr = isSynced(tpr);
    close(tpr);
elseif strcmp(scpr,'Imager')
    [name, date] = fetch1( Scans(key),'file_name','exp_date' );
    if ~isempty(strfind(name,'.h5'));fend = [];else fend = '.h5';end
    filename = getLocalPath(['M:/IntrinsicImaging/' datestr(date, 'yymmdd') '/' name fend]);
    isSyncedTpr = 0;
else
    isSyncedTpr = 0;
end

override = 0;
% isSyncedTpr = 0;
% synchronize all the scans from that day
% isempty(VisStims(rmfield(key,'scan_idx')))
if  ~isSyncedTpr && state && strcmp(scpr,'MPScan') && ...
        ~isempty(Scans(key,'exp_date > "2012-05-01"')) && ~override 
    syncStim2TpDay(key.exp_date,'stimPath',stimPath)
    isSyncedTpr = 1;
end
% isSyncedTpr = 0;
if state && strcmp(scpr,'MPScan') && ~override
    stimfilename = syncStim2Tpr(filename(1:end-6),'stimPath',stimPath);
    tpr = tpReader( Scans(key) );
    if isSynced(tpr)
        tuple.frame_timestamps = readSyncedTimestamps(tpr);
    else
        return
    end
elseif  state && strcmp(scpr,'MPScan') && override %override && state && strcmp(scpr,'MPScan')
    stimfilename = syncStim2Tpr(filename(1:end-6),'stimPath',stimPath,'manual',1);
    tpr = tpReader( Scans(key) );
    if isSynced(tpr)
        tuple.frame_timestamps = readSyncedTimestamps(tpr);
    else
        return
    end
elseif ~isSyncedTpr && ~state && strcmp(scpr,'MPScan')
    tpr = tpReader( Scans(key) );
    try
        synchronize(tpr);
    catch
        return
    end
    stimfilename = findStimFiles(filename(1:end-6),'stimPath',stimPath);
    tpr = tpReader( Scans(key) );
    tuple.frame_timestamps = readSyncedTimestamps(tpr);
elseif strcmp(scpr,'AOD') && state
    try
        fps = fetch1(Movies(key),'fps');
    catch
        fps = 40;
    end
    [stimfilename, tuple.frame_timestamps] = syncStim2AOD(filename, ...
        'stimPath',stimPath,'scprog',scpr,'date',fetch1(Scans(key),'exp_date'),...
        'fps',fps,'aodReader',aodReader,'manual',1);
elseif strcmp(scpr,'Unirec') && state
    [stimfilename, tuple.frame_timestamps] = syncStim2Unirec(filename, ...
        'stimPath',stimPath,'scprog',scpr,'date',fetch1(Scans(key),'exp_date'));
elseif strcmp(scpr,'Imager') && state
    [stimfilename, tuple.frame_timestamps] = syncStim2Imager(filename, ...
        'stimPath',stimPath,'scprog',scpr,'date',fetch1(Scans(key),'exp_date'),'manual',1);
elseif strcmp(scpr,'Imager') && ~state
    stimfilename = {'other'};
    tuple.frame_timestamps = [];
elseif strcmp(scpr,'ScanImage') && state
    tpr = tpReader( Scans(key) );
    [stimfilename, tuple.frame_timestamps] = syncStim2ScanImage(tpr);
else
    stimfilename = findStimFiles(filename(1:end-6),'stimPath',stimPath);
    tpr = tpReader( Scans(key) );
    tuple.frame_timestamps = readSyncedTimestamps(tpr);
end

if isempty( stimfilename )
    disp('No matching stimulus file is found' );
else
    for iStim = 1:length(stimfilename)
        tuple.stim_idx = iStim;
        tuple.stim_filename = stimfilename{iStim};
        [foo, tuple.exp_type] = fileparts( tuple.stim_filename );  %#ok<ASGLU> % the *.mat file must be named after the experiment type
        if ~strcmp(tuple.stim_filename,'other')
            tmp = load(getLocalPath( tuple.stim_filename ));
            stim = tmp.stim;
            stim.params.trials = rmfield(stim.params.trials,'swapTimes');
            tuple.stim_file = stim;
        else
            tuple.stim_file = [];
        end
        % insert things
        insert( obj, tuple );
    end
end
