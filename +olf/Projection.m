%{
-> reso.FluorescenceTrace
---
projection                    : binary                    # mask with signal in red channel
%}


classdef Projection < dj.Computed
    
    properties
         keySource = (reso.FluorescenceTrace & olf.Sync)
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
            
            % get avg image
            im = fetch1(reso.SummaryImagesAverage & key,'average_image');

            % get mask
            [px,wt] = fetch1(reso.SegmentationMask & key,'pixels','weights');

            BW = backim;
            BW(px) = wt;

            B = bwboundaries(BW,'noholes');
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