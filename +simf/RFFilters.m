%{
->simf.RFParams
filter_id                             : smallint      #
---
filter                            : blob                          #
%}

classdef RFFilters < dj.Manual
    methods
        function IM = plot(self)
            [filters] = fetchn(self,'filter');
            filters = cell2mat(cellfun(@(x) x(:),filters,'uni',0)');
            im = reshape(permute(reshape(filters,16,16,16,16),[1 3 2 4]),256,[]);
            if nargout
                IM = im;
            else
                figure
                imagesc(imresize(im,10))
            end
            axis image
            axis off
            colormap(cbrewer('div','RdBu',100))
        end
        
         function IM = plotSpace(self,upsample)
             
             if nargin<2
                 upsample = 1;
             end
             
            [filters] = fetchn(self,'filter');
            filts = cell2mat(cellfun(@(x) x(:),filters,'uni',0)');
            mx = max(filts(:));
            mn = min(filts(:));
            figure
            for ifilt = 1:length(filters)
                subplot_tight(sqrt(length(filters)),sqrt(length(filters)),ifilt,0.002)
                 imagesc(imresize(filters{ifilt},upsample),[mn,mx])
                 set(gca,'xtick',[],'ytick',[],'XColor',[0.85 0.85 0.85],'ycolor',[0.85 0.85 0.85])
                 axis square
            end
            colormap(cbrewer('div','RdBu',100))
            set(gcf,'position',[200 200 500 500])
            shg
        end
    end
end

