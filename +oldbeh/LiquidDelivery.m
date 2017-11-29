%{
beh.LiquidDelivery (manual) # 
-> beh.Session
timestamp       : bigint                 # 
---
pulse_time=null             : smallint                      # 
microl_per_pulse=null       : float(6,2)                    # 
liquid_type=null            : enum('apple_juice','water')   # 
%}


classdef LiquidDelivery < dj.Relvar
        methods
        function self = Session(varargin)
            self.restrict(varargin{:})
        end
        
        function plotDelivery(self)
            figure

%                 mice = fetch(vis2p.Mice & (beh.LiquidDelivery & ['timestamp>' num2str(timestamp(floor(now)))]));
            mice = fetch(vis2p.Mice & self);

            for imouse = 1:length(mice);
               
                [times, pulse] = fetchn(self & mice(imouse),'timestamp','microl_per_pulse');
                [days, ~, IC]= unique(floor(timestamp(times,'matlab')),'rows');
                water = nan(length(days),1);
                
                for iday = 1:length(days);
                   water(iday) = roundall(sum(pulse(IC==iday))/1000,0.1);
                    
                end
                subplot(ceil(sqrt(length(mice))),floor(sqrt(length(mice))),imouse)
                plot(water);
%                 set(gca,'xtick',1:length(days),'xticklabel',datestr(days),...
%                     'XTickLabelRotation',90,'box','off')
                   set(gca,'xtick',1:length(days),'xticklabel',datestr(days),...
                    'box','off')
                hold on
                plot([0 length(days)+1],[1 1],'r-.')
                xlim([0 length(days)])
                ylim([0 3])
                title(num2str(mice(imouse).mouse_id))
            end
        end
    end

end