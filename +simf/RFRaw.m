%{
->simf.RFParams
->stimulus.Movie
---
%}


classdef RFRaw < dj.Computed
    properties
        keySource  = (simf.RFParams & simf.RFFilters) * stimulus.Movie
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            target_fps = 2; % ~500ms bins
            [x_sz, y_sz] = fetch1(simf.RFParams & key,'x_size','y_size');
            [filters, keys] = fetchn(simf.RFFilters & key,'filter');
            filters = cell2mat(cellfun(@(x) x(:),filters,'uni',0)');
            movie_keys = fetch(stimulus.MovieClip & key);
            responses = cell(length(movie_keys),1);
            parfor imovie=1:length(movie_keys)
                
                filename = export(stimulus.MovieClip & movie_keys(imovie), 'temp');
                
                vr = VideoReader(filename{1});
               
                frame_rate = vr.FrameRate;
                skip_frames = round(frame_rate/target_fps);
                file = nan(vr.Height,vr.Width,3,floor(vr.Duration*frame_rate));
                for i = 1:size(file,4)
                    file(:,:,:,i) = vr.readFrame;
                end
                file = file(:,:,:,1:skip_frames:end)/255;
                delete(filename{1})
                
                file = imresize(squeeze(file(:,:,1,:)),y_sz/vr.Height);
                file = file(1:y_sz,round(size(file,2)/2 - x_sz/2)+1:round(size(file,2)/2+ x_sz/2),:);
                file = reshape(file,size(filters,1),size(file,3));
                
%                 responses(:,end+1:end+size(file,2)) = filters'*file;
                responses{imovie} = filters'*file;
            end
            responses = cell2mat(responses');
            insert(self,key)
           [keys.movie_name] = deal(key.movie_name);
            for ikey = 1:length(keys)
                rf_key = keys(ikey);
                rf_key.response = responses(ikey,:);
                makeTuples(simf.RFRawResponses,rf_key)
            end
        end
    end
end

