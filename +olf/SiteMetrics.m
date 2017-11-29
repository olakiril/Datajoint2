%{
olf.SiteMetrics (computed) #
-> olf.RespOpt
-> preprocess.Spikes
-> preprocess.MaskType
---
sparseness_on              : float                    # sparseness
sparseness_off             : float                    # sparseness
pactive_on                 : float                    # ratio of active cells
pactive_off                : float                    # ratio of active cells
signal_corr                : mediumblob               # signal correlations
noise_corr                 : mediumblob               # noise correlations
total_corr                 : mediumblob               # total correlations
eu_dist                    : mediumblob               # distances between cells
%}


classdef SiteMetrics < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel =  (preprocess.Spikes * preprocess.MaskType  & preprocess.MaskClassificationMaskType) * olf.RespOpt
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
 
            
            [resp_on, resp_off] = fetchn(olf.OlfResponses * preprocess.MaskClassificationMaskType & key,'resp_on','resp_off');
            resp_on = cell2mat(cellfun(@(x) permute(x,[3 1 2]),resp_on,'uni',0));
            resp_off = cell2mat(cellfun(@(x) permute(x,[3 1 2]),resp_off,'uni',0));
            tuple = key;
            
            if isempty(resp_on); display('No traces!');return; end
            
            % sparseness
            tuple.sparseness_on = nanmean(sparseness(resp_on(:,:)));
            tuple.sparseness_off = nanmean(sparseness(resp_off(:,:)));
            
            % p active
            mn = repmat(nanmean(resp_on(:,:),2),1,size(resp_on,2),size(resp_on,3));
            sd = repmat(nanstd(resp_on(:,:),[],2),1,size(resp_on,2),size(resp_on,3));
            tuple.pactive_on = mean(resp_on(:)>mn(:)+sd(:));
            mn = repmat(nanmean(resp_off(:,:),2),1,size(resp_off,2),size(resp_off,3));
            sd = repmat(nanstd(resp_off(:,:),[],2),1,size(resp_off,2),size(resp_off,3));
            tuple.pactive_off = mean(resp_off(:)>mn(:)+sd(:));

            % correlations
            traces = fetchn(preprocess.SpikesRateTrace & key & (preprocess.MaskClassificationMaskType & key),'rate_trace');
            traces = (cell2mat(traces'));
            tuple.total_corr = corr(traces);
            tuple.signal_corr = corr(nanmean(resp_on,3)');
            z = zscore(resp_on,[],3);
            tuple.noise_corr = corr(z(:,:)');
            
            % distances
            [x,y,z] = fetchn(preprocess.MaskCoordinates * preprocess.MaskClassificationMaskType & key,...
                'xloc','yloc','zloc');
            tuple.eu_dist = squareform(pdist([x,y,z]));
            
            % insert
            insert( obj, tuple );
            
        end
    end
end