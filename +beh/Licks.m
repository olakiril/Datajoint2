%{
beh.Licks (manual) #
-> beh.Session
timestamp       : bigint                 #
---
%}


classdef Licks < dj.Relvar
    
    
    methods
        function self = Traces(varargin)
            self.restrict(varargin{:})
        end
        
        function plot(self,tmst)
            
            figure
            set(gcf,'name','Correct Stimuli')
            k = [];
            if nargin<2
                tmst = 'session_timestamp>"2016-07-15 09:00:00"';
            end
%             tmst = opts;
            mice = unique(fetchn(beh.Session & tmst & self & 'exp_type > "Freerun"','mouse_id'));
            for ii = 1:length(mice)
                
                k.mouse_id = mice(ii);
                p_names = unique(fetchn(beh.Session & k & tmst & 'exp_type>"Freerun"','stimuli'));
                p_types = unique(fetchn(beh.Session & k & tmst & 'exp_type>"Freerun"','rewarded_stimuli'));
                k.period_type = p_types{end};
                
                wtimes = double(fetchn(beh.StimPeriods & k & tmst,'timestamp'));
                ltimes = double(fetchn(beh.Licks & k & tmst,'timestamp'));
                stime = min([wtimes(:);ltimes(:)]);
                
                wtimes = (wtimes-stime)/1000/60;
                ltimes = (ltimes-stime)/1000/60;
                subplot(1,length(mice),ii)
                for i = 1:length(wtimes);
                    if i ~=length(wtimes)
                        licks = ltimes(ltimes>wtimes(i) & ltimes<wtimes(i+1));
                    else
                        licks = ltimes(ltimes>wtimes(i));
                    end
                    
                    plot(licks-wtimes(i),ones(size(licks))*i,'.k')
                    hold on
                end
                set(gca,'ydir','reverse','box','off')
                xl = get(gca,'xlim');
                xlim([-1 5])
                ylim([0 length(wtimes)+1])
                plot([0 0],[1 length(wtimes)],'b')
                plot([1 1],[1 length(wtimes)],'--r')
                if ii==1
                    ylabel('Trial #')
                end
                xlabel('Time (min)')
                title([num2str(mice(ii))])
                set(gca,'ytick',[])
            end
            
            figure
            set(gcf,'name','Incorrect Stimuli')
            k = [];
            for ii = 1:length(mice)
                
                k.mouse_id = mice(ii);
                
                p_names = unique(fetchn(beh.Session & k & tmst & 'exp_type>"Freerun"','stimuli'));
                periods =  strsplit(p_names{end},',');
                p_types = unique(fetchn(beh.Session & k & tmst & 'exp_type>"Freerun"','rewarded_stimuli'));
                k.period_type = p_types{end};
                
                % wrong stimulus
                idx = ~strcmp(periods,k.period_type);
                if sum(idx)~=0
                    k.period_type = periods{idx};
                else
                    continue
                end
                
                wtimes = double(fetchn(beh.StimPeriods & k & tmst,'timestamp'));
                ltimes = double(fetchn(beh.Licks & k & tmst,'timestamp'));
                stime = min([wtimes(:);ltimes(:)]);
                
                wtimes = (wtimes-stime)/1000/60;
                ltimes = (ltimes-stime)/1000/60;
                subplot(1,length(mice),ii)
                for i = 1:length(wtimes);
                    if i ~=length(wtimes)
                        licks = ltimes(ltimes>wtimes(i) & ltimes<wtimes(i+1));
                    else
                        licks = ltimes(ltimes>wtimes(i));
                    end
                    
                    plot(licks-wtimes(i),ones(size(licks))*i,'.k')
                    hold on
                end
                set(gca,'ydir','reverse','box','off')
                xl = get(gca,'xlim');
                xlim([-1 5])
                ylim([0 length(wtimes)+1])
                plot([0 0],[1 length(wtimes)],'b')
                plot([1 1],[1 length(wtimes)],'--r')
                if ii==1
                    ylabel('Trial #')
                end
                xlabel('Time (min)')
                title([num2str(mice(ii))])
                set(gca,'ytick',[])
            end
            
        end
        
        function plotLinear(self,k)
            licks = fetchn(beh.Licks & k,'timestamp');
            [periods,periodst]=fetchn(beh.StimPeriods & k,'period_type','timestamp');
            water = fetchn(beh.LiquidDelivery & k,'timestamp');
            air = fetchn(beh.AirDelivery & k,'timestamp');

            plot((periodst(strcmp(periods,'obj1_default0004.png'))-periodst(1))/1000,1,'.g')
            hold on;plot((periodst(strcmp(periods,'obj2_default0002.png'))-periodst(1))/1000,1,'.r')
            hold on;plot((periodst(strcmp(periods,'endTrial'))-periodst(1))/1000,1,'.k')
            hold on;plot((periodst(strcmp(periods,'punishment'))-periodst(1))/1000,2,'.r')
            hold on;plot((water-periodst(1))/1000,1.5,'.g')
            hold on;plot((licks-periodst(1))/1000,0,'.k')

            ylim([-1 10])
            grid on
        end
    end
    
    
end