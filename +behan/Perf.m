%{
# Performance with 2AFC tasks
-> behan.SesPerf
trial_idx       : smallint               # unique condition index
---
p                    : int                            # hit trials
%}

classdef Perf < dj.Computed
    properties
        keySource = beh.Trial & behan.SesPerf
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
        function  plot(self)
            %             mice = fetchn(beh.SetupInfo & self & 'animal_id>0','animal_id');
            
            figure; s= [];
            [Perf_out, Days_out, mice, Trials_out] = getPerformance(self);
            for imouse = 1:length(mice)
                s(imouse) = subplot(round(sqrt(length(mice))),ceil(sqrt(length(mice))),imouse);
                yyaxis left
                plot([1:length(Days_out{imouse})],Perf_out{imouse},'linewidth',2)
                hold on
                plot([1 length(Days_out{imouse})],[0.5 0.5],'-.','color',[0.6 0.6 0.6],'linewidth',2)
                set(gca,'xtick',1:length(Days_out{imouse}),'xticklabel',datestr(Days_out{imouse},'mm/dd'),'XTickLabelRotation',45)
                ylabel('Performance')
                yyaxis right
                plot([1:length(Days_out{imouse})],Trials_out{imouse})
                ylabel('Trials')
                title(sprintf('%d',mice(imouse)))
                xlim([1 length(Days_out{imouse})])
                grid on
                
            end
            
            linkaxes(s,'y')
        end
        
        function [Perf_out, Days_out, mice, Trials_out] = getPerformance(self)
            mice = unique(fetchn(self & 'animal_id>0','animal_id'));
            
            Days_out = []; Perf_out = [];Trials_out = [];
            for imouse = 1:length(mice)
                try

                    sess_keys = fetch(beh.Session & self & 'exp_type="CenterPort" OR exp_type="CenterPortTrain"' & sprintf('animal_id = %d',mice(imouse)));
                    
                    perf = [];times = [];
                    for isession = 1:length(sess_keys)
                        [perf{isession},t] = fetchn(behan.Perf * beh.Trial & sess_keys(isession),'p','start_time');
                        times{isession} =  msec2tmst(beh.Session & sess_keys(isession),t);
                    end
                    perf = cell2mat(perf');
                    times = cell2mat(times');
                    
                    days = unique(floor(times));
                    Perf = nan(length(days),1);trial_num = Perf;
                    for iday = 1:length(days)
                        idx = times>days(iday) & times<days(iday)+1;
                        Perf(iday) = nanmean(perf(idx));
                        trial_num(iday) = sum(idx);
                    end
                    
                    Perf_out{imouse} = Perf;
                    Days_out{imouse} = days;
                    Trials_out{imouse} = trial_num;
                    
                end
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