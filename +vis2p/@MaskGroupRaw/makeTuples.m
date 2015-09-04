function makeTuples( obj, key )

import vis2p.*

[path,name,scan_prog] = fetch1( Experiments*Scans(key), 'directory','file_name','scan_prog');

if strcmp(scan_prog,'AOD')
    if ~isempty(Scans(key,'exp_date > "2012-05-01"'))
        dirs = dir(['M:\Mouse\' name(1:10) '*']);
        names = vertcat(dirs.name);
        timediff = str2num(names(:,[12 13 15 16 18 19]))- str2num( name([12 13 15 16 18 19])); %#ok<ST2NM>
        timeind = find(timediff<0);
        [~,i] = max(timediff(timeind));
        itime = timeind(i);
        filename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' name];
        [traces,timestamps, coordinates, traces2,fps] = loadAODTraces(filename,0);
    end
    
    otherTimes = fetch1(VisStims(key,'stim_idx = 1'),'frame_timestamps');
    timestamps = (timestamps - timestamps(1))*1000+otherTimes(1);
    masknums = fetchn(MaskTraces(key),'masknum');
else
    timestamps = fetch1(VisStims(key,'stim_idx = 1'),'frame_timestamps');
    fps = fetch1(Movies(key),'fps');
    [traces,traces2,masknums] = fetchn(MaskTraces(key),'calcium_trace','red_trace','masknum');
    [coordinates(:,1), coordinates(:,2), coordinates(:,3)] = fetchn(MaskCells(key),'img_x','img_y','img_z');
    traces = cat(2,traces{:});
    traces2 = cat(2,traces2{:});
end

tuple = key;
tuple.frame_timestamps = timestamps;
tuple.fps = fps;

insert( obj, tuple );
makeTuples(vis2p.MaskTracesRaw,key,traces,traces2,coordinates,masknums)

