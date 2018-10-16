%{
rf_opt       : smallint unsigned      #
---
x_size                            : smallint                            # x filter size
y_size                            : smallint                            # y filter size
type                              : enum('pixels','ICA','biICA','multiICA')        # filter type
resize_factor                     : double
%}


classdef RFParams < dj.Lookup
    methods
        function createFilters(self)
            [type,x_sz,y_sz,rsz] = fetch1(self,'type','x_size','y_size','resize_factor');
            switch type
                case 'pixels'
                    filters = reshape(diag(ones(x_sz*y_sz,1)),y_sz,x_sz,x_sz*y_sz);
                case 'ICA'
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                case 'biICA'
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                    c = corr(abs(filt.A));
                    c(logical(diag(ones(size(c,1),1))))=nan;
                    for ifilter = 1:size(filters,3)
                        [~,idx] = nanmax(c(ifilter,:));
                        filters(:,:,ifilter) = filters(:,:,ifilter) + filters(:,:,idx);
                    end
                case 'mmICA'
                    [responses, ids] = fetchn(simf.RFRawResponses & 'movie_name = "MadMax"' & 'rf_opt = 2','response','filter_id');
                    [~,idx]= sort(ids);
                    responses = cell2mat(responses(idx));
                    c = corr(responses');
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                    for ifilter = 1:size(filters,3)
                        [~,idx] = sort(c(ifilter,:),'descend');
                        filters(:,:,ifilter) = mean(filters(:,:,idx(1:2)),3);
                        c(idx(2),idx(1)) = 0;
                    end
                case 'obj1ICA'
                    [responses, ids] = fetchn(simf.RFRawResponses & 'movie_name = "obj1v6"' & 'rf_opt = 2','response','filter_id');
                    [~,idx]= sort(ids);
                    responses = cell2mat(responses(idx));
                    c = corr(responses');
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                    for ifilter = 1:size(filters,3)
                        [~,idx] = sort(c(ifilter,:),'descend');
                        filters(:,:,ifilter) = mean(filters(:,:,idx(1:2)),3);
                        c(idx(2),idx(1)) = 0;
                    end
                case 'objICA'
                    Resp = [];
                    [responses, ids] = fetchn(simf.RFRawResponses & 'movie_name = "obj1v6"' & 'rf_opt = 2','response','filter_id');
                    [~,idx]= sort(ids);
                    Resp{1} = cell2mat(responses(idx));
                    [responses, ids] = fetchn(simf.RFRawResponses & 'movie_name = "obj2v6"' & 'rf_opt = 2','response','filter_id');
                    [~,idx]= sort(ids);
                    Resp{2} = cell2mat(responses(idx));
                    [responses, ids] = fetchn(simf.RFRawResponses & 'movie_name = "obj3v6"' & 'rf_opt = 2','response','filter_id');
                    [~,idx]= sort(ids);
                    Resp{3} = cell2mat(responses(idx));
                    [responses, ids] = fetchn(simf.RFRawResponses & 'movie_name = "obj4v6"' & 'rf_opt = 2','response','filter_id');
                    [~,idx]= sort(ids);
                    Resp{4} = cell2mat(responses(idx));
                    responses = cell2mat(Resp);
                    c = corr(responses');
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                    for ifilter = 1:size(filters,3)
                        [~,idx] = sort(c(ifilter,:),'descend');
                        filters(:,:,ifilter) = mean(filters(:,:,idx([1 5])),3);
                        c(idx(2),idx(1)) = 0;
                    end
            end
            
            if rsz>1
                for ifilter = 1:size(filters,3)
                    filt = imresize(filters(:,:,ifilter),rsz);
                    im = abs(filters(:,:,ifilter));
                    im2 = imresize(normalize(imfilter(im,gausswin(7)*gausswin(7)','same')),2);
                    Ibw = imfill(imbinarize(im2),'holes');
                    stat = regionprops(logical(Ibw),'centroid');
                    x_off = round((size(im,1)/2 - stat(1).Centroid(1)/2)/size(im,1)*(size(filt,1) - size(im,1)));
                    y_off = round((size(im,2)/2 - stat(1).Centroid(2)/2)/size(im,2)*(size(filt,2) - size(im,2)));
                    filters(:,:,ifilter) = filt((1:size(filt,1)/rsz)+(size(filt,1) - size(im,1))/2-y_off,(1:size(filt,2)/rsz)+(size(filt,2) - size(im,2))/2-x_off);
                end
            end
            
            % insert filters
            key = fetch(self);
            for ifilter = 1:size(filters,3)
                key.filter_id = ifilter;
                key.filter = filters(:,:,ifilter);
                insert(simf.RFFilters, key)
            end
        end
     
    end
end

