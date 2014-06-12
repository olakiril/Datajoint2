function makeTuples(obj,key)

import vis2p.*

if strcmp(fetch1(Scans(key),'scan_prog'),'Unirec')
    if key.trace_opt == 6
    [path,name] = fetch1(Experiments*Scans(key), 'directory','file_name' );
    filename = getLocalPath([path '/' name '%d.h5']);
     br = baseReader(filename);
     
     % take the 2norm and bin;
    filter = filterFactory.createBandpass(400, 600, 5800, 6000, getSamplingRate(br));
    fr = filteredReader(br, filter);
    traces = sum(fr(:,1:4).^2,2);
    traces = sqrt(traces);

    % crude filter at 10 hz
    d = max(1,round(0.1*getSamplingRate(br)));
    k = ones(d,1)/d;
    traces = convn(traces,k,'same');
    
    tuple = key;
    tuple.masknum = 1;
    tuple.trace = single(traces);
    tuple.skewness = skewness( traces );
    insert( obj, tuple );
    else
        return
    end
else
    % get traces and frames
    if strcmp(fetch1(Scans(key),'scan_prog'),'AOD')
        [masknums,traces(:,:,1),traces(:,:,2)] = ...
            fetchn(MaskTraces(key,'masknum <> -1 and masknum <> -3 and masknum <> 0'),...
            'masknum','calcium_trace','red_trace');
    elseif  strcmp(fetch1(Scans(key),'scan_prog'),'ScanImage')
          [masknums,traces] = ...
            fetchn(MaskTraces(key,'masknum <> -1 and masknum <> -3 and masknum <> 0'),...
            'masknum','calcium_trace');
    else
        [masknums,traces(:,:,1),traces(:,:,2),traces(:,:,3)] = ...
            fetchn(MaskTraces(key,'masknum <> -1 and masknum <> -3 and masknum <> 0'),...
            'masknum','calcium_trace','red_trace','annulus_trace');
    end
    fps = fetch1(Movies(key),'fps');  %#ok<NASGU>
    traceOpt = fetch(TracesOpt(key),'*'); %#ok<NASGU>
    options = fetch1(TracesOpt(key),'trace_computation');
    options = strread(options, '%s','delimiter',',');
    traceOpt.key = key;
    
    % filter traces
    traceTypes = size(traces,3);
    traces = double([traces{:}]);  % traces are in columns
    traces = reshape(traces,size(traces,1),[],traceTypes);
    for iopt = 1:length(options)
        [traces, qual, traceOpt]= eval([options{iopt} '(traces,fps,traceOpt)']);
    end
    
    % insert cell traces and annulus traces
    for imask=1:length(masknums)
        tuple = key;
        tuple.masknum = masknums(imask);
        tuple.trace = single(traces(:,imask));
        tuple.quality = qual(imask);
        insert( obj, tuple );
    end     
    
    if  strcmp(fetch1(Scans(key),'aim'),'patching')&&  ~isempty(EphysTraces(key))
        tuple = key;
        tuple.masknum = 0;
 
        % construct ephys trace with frame times
        imsize = fetch1(Movies(key),'nframes');
        spikeTimes = fetch1(EphysTraces(key),'spike_times');
        eTrace = zeros(imsize,1);
        bTrace = ceil(spikeTimes*fps/1000);
        bTrace(bTrace>imsize) = imsize;
        for iSpike = 1:length(bTrace)
            eTrace(bTrace(iSpike)) = eTrace(bTrace(iSpike)) + 1;
        end
        
        tuple.trace = single(eTrace);
        insert( obj, tuple );
    end
end


