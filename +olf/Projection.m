%{
-> reso.FluorescenceTrace
---
projection                    : int8                    # mask with signal in red channel
%}


classdef Projection < dj.Computed
    
    properties
         keySource = reso.FluorescenceTrace & (reso.MaskClassificationType & 'type = "tufted"') & olf.Sync
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
            
            % get avg image
            key2 = key;
            key2.channel = 2;
            im = fetch1(reso.SummaryImagesAverage & key2,'average_image');
            backim = false(size(im));
            
            % get mask
            [px,wt] = fetch1(reso.SegmentationMask & key,'pixels','weights');

            BW = backim;
            BW(px) = wt;
              
            stat = regionprops(BW,'EquivDiameter');
            BW2 = imerode(BW,strel('sphere',round(stat.EquivDiameter/3)));
            BWD = convn(BW,gausswin(40)*gausswin(40)','same');
            mask_sig= ttest2(im(BW2),im(BWD & ~BW),'tail','right','Alpha',0.001);

            % insert
            tuple = key;
            tuple.projection = mask_sig;
            insert( obj, tuple );
        end
    end
end