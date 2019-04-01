%{
# Performance with 2AFC tasks
animal_id       : int                    # animal id
session_id      : smallint               # session number
trial_idx       : smallint               # unique condition index
---
p                    : int                            # hit trials
%}

classdef Perf < dj.Computed
    properties
        keySource = proj(beh.Trial) & proj(beh.Session & (beh.RewardCond & 'probe=1') & (beh.RewardCond & 'probe=2') & 'animal_id>0')
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % fetch session time
            [start, stop, probe] = fetch1(beh.Trial * beh.RewardCond & key,'start_time','end_time','probe');
            
%             % fetch session time
%             [start, stop, probe] = fetch1(beh.Trial * beh.MovieClipCond * beh.RewardCond & key,'start_time','end_time','probe');
            
            % get licks
            resp_probe = fetchn(beh.Lick & key & sprintf('time>%d',start) & sprintf('time<%d',stop),'probe','ORDER BY time');
            
            
            tuple.animal_id = key.animal_id;
            tuple.session_id = key.session_id;
            tuple.trial_idx = key.trial_idx;
            if isempty(resp_probe) && probe==0
                tuple.p = 1;
            elseif isempty(resp_probe) && probe~=0
                tuple.p = nan;
            elseif resp_probe(1)==probe
                tuple.p = 1;
            else
                tuple.p = 0;
            end
            
            % insert
            insert(self,tuple)
            
        end
    end
    
    methods
        function plot(self)
            %             mice = fetchn(beh.SetupInfo & self & 'animal_id>0','animal_id');
            mice = unique(fetchn(self & 'animal_id>0','animal_id'));
            
            figure
            for imouse = 1:length(mice)
                subplot(ceil(sqrt(length(mice))),ceil(sqrt(length(mice))),imouse)
                [sessions,session_tmst] = fetchn(beh.Session & 'exp_type="CenterPort"' & sprintf('animal_id = %d',mice(imouse)),'session_id','session_tmst');
                perf = [];
                for isession = 1:length(sessions)
                    perf(isession) = nanmean(fetchn(behan.Perf  & sprintf('session_id = %d',sessions(isession))  & sprintf('animal_id = %d',mice(imouse)),'p')) ;
                end
                plot(perf)
                hold on
                plot([1 isession],[0.5 0.5],'-.','color',[0.6 0.6 0.6])
                set(gca,'xtick',1:isession,'xticklabel',datestr(session_tmst),'XTickLabelRotation',45)
                title(sprintf('%d',mice(imouse)))
                xlim([1 isession])
                ylim([0.4 0.8])
            end
        end
        
        function perf = getSessPerformance(self)
            sessions = fetch(beh.Session & self);
            perf = [];
            for isession = 1:length(sessions)
                perf(isession) = nanmean(fetchn(behan.Perf  & sessions(isession),'p')) ;
            end
        end
    end
end