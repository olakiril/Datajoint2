%{
#
-> experiment.Scan
---
fps               : double               # frames per second of the recording
frames            : double               # total frames recorded
trace_fs          : double               # final sampling rate of the respiration signal 
trace             : mediumblob           # raw respiration signal
%}


classdef RespRaw < dj.Imported
    
    properties
        keySource = experiment.Scan & experiment.AutoProcessing & proj(olf.Session & 'session_timestamp>"2018-08-15"','mouse_id->animal_id')
    end
    
    methods (Access=protected)
        function makeTuples(self, key)
                        
            [file,file_name] = fetch1(olf.Session & key,'file_name','file_name');
            file_key =[];
            file_key.animal_id = key.mouse_id;
            file_key.filename = file_name;
            path = fetch1(experiment.Session & (experiment.Scan & file_key),'scan_path');
            filetype = getLocalPath(fullfile(path,sprintf('%s%s',file,'*.tif')));
            display(['Reading file: ' filetype])
            k = [];
            k.animal_id = key.mouse_id;
            k.filename = file;
            [key.session, key.scan_idx, key.animal_id] = ...
                fetch1(experiment.Scan & k,'session','scan_idx','animal_id');
            
            % get data in chunks
            reader = ne7.scanimage.Reader5(filetype);
            sz = size(reader);
            frame_chuck = 1000;
            data_phd = nan(sz(1),sz(4),sz(5));
            for i = 1:ceil(sz(5)/frame_chuck)
                idx = 1:frame_chuck + (i-1)*frame_chuck;
                idx(idx>sz(5)) = [];
                data_phd(:,:,idx) = squeeze(nanmean(reader(:,:,find(reader.channels==4),:,idx),2));
            end
            
            % get the sampling close to 100Hz
            y_pixel_bin = round(round(100/(reader.fps * sz(4)))/2)*2;
            assert(mod(sz(1),y_pixel_bin)==0,sprintf('%d pixels not devided by %d',sz(1),y_pixel_bin))
            key.trace = reshape(nanmean(reshape(data_phd,sz(1)/y_pixel_bin,y_pixel_bin,[]),1),[],1);
            key.trace_fs = reader.fps*sz(4)*4;
            key.fps = reader.fps;
            key.frames = sz(5);
            
            % insert 
            insert(self, key);
        end
    end

end