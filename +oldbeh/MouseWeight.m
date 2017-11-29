%{
beh.MouseWeight (manual) #
-> vis2p.Mice
timestamp       : timestamp              #
---
weight                      : double(5,2)                   #
%}


classdef MouseWeight < dj.Relvar
    
    methods
        function self = MouseWeight(varargin)
            self.restrict(varargin{:})
        end
        
        function plot(self)
            mice = unique(fetchn(self,'mouse_id'));
            figure
            set(gcf,'name','Mouse Weight','NumberTitle','off',...
                'units','normalized','outerposition',[0 0 1 1])

            k = [];
            for imouse = 1:length(mice)
                k.mouse_id = mice(imouse);
                [weights, times] = fetchn(self & k,'weight','timestamp');
                
                subplot(round(sqrt(length(mice))),ceil(sqrt(length(mice))),imouse)
                plot((weights-weights(1))*100/weights(1))
                set(gca,'xtick',1:length(weights),'xticklabel',...
                    cellfun(@(x) x(1:10),times,'uni',0),'xticklabelrotation',45);
                hold on
                plot([0 length(weights)+1],[-20 -20],'-.r')
                
                title(num2str(mice(imouse)))

                if imouse==1;ylabel('% weight change');end
                
            end
        end
    end
end