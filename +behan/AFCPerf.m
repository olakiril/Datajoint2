%{
# Performance with 2AFC tasks for each day
->beh.Session
day                   : timestamp                      # day starts at 12AM
---
ht                    : int                            # hit trials
fa                    : int                            # false alarm trials
ms                    : int                            # miss trials
cr                    : int                            # correct rejection trials
pval                  : float                          # bootstrap p value of HIT_RATE > FA_RATE
dprime                : float                          # Sensitivity d? = z(Hit) ? z(FA)
perf                  : float                          # performance (HITS + CR/ALL TRIALS)
dconf                 : float                          # 95% d' Confidence
%}

classdef AFCPerf < dj.Computed
    properties
        keySource = beh.Session & (beh.RewardCond & 'probe=1') & (beh.RewardCond & 'probe=2') & (beh.MovieClipCond & ...
            'movie_name="obj1v4" OR movie_name="obj2v4" OR movie_name="obj3v4" OR movie_name="obj3v4"') & 'animal_id>0'
    end
    
    methods
        
        function plot(self, type)
            
            figure
            if nargin<2; type = 'perf';end
            
            allmice = fetch(mice.Mice & self);
            sh = nan(length(allmice),1);
            for imouse =1:length(allmice)
                mouse = allmice(imouse);
                sh(imouse) = subplot(ceil(sqrt(length(allmice))),ceil(sqrt(length(allmice))),imouse);
                switch type
                    case 'perf'
                        [ht,fa,ms,cr,days] = fetchn(self & mouse,'ht','fa','ms','cr','day');
                        day = unique(days);
                        val = nan(length(day),1);
                        for iday = 1:length(day)
                            idx = strcmp(days,day{iday});
                            val(iday) = (sum(ht(idx))+sum(cr(idx)))/sum(sum(ht(idx))+sum(fa(idx))+sum(ms(idx))+sum(cr(idx)));
                        end
                        line = 0.5;
                    case 'dprime'
                        line = 0;
                        [val,day] = fetchn(self & mouse,'dprime','day');
                end
                
                dn = datetime(day);
                plot(dn,val)
                hold on
                xlim([min(dn) max(dn)])
                plot([min(dn) max(dn)],[1 1]*line,'-.','color',[0.5 0.5 0.5])
                title(sprintf('Mouse ID: %d',mouse.animal_id))
                if imouse==length(allmice)
                    %                     set(gca,'xtick',min(dn):2:max(dn),'XTickLabel',datestr(datevec(min(dn):2:max(dn)),'DD-mm'),'XTickLabelRotation',45)
                    %set(gca,'xtick',min(dn):max(dn))
                    %datetick('x','dd-mm','keepticks','keeplimits')
                    xtickformat('dd-MM')
                    set(gca,'xticklabelrotation',45)
                else
                    set(gca,'xticklabel',[])
                end
            end
            linkaxes(sh,'x')
        end
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % fetch session times
            [start_times, session_times, end_times, rprobe] = ...
                fetchn((beh.Trial & key) * beh.Session * (beh.Condition & key) * (beh.RewardCond & key),...
                'start_time','session_tmst','end_time','probe','ORDER BY start_time');
            session_times = datenum(session_times,'YYYY-mm-dd HH:MM:SS');
            
            % start day at 12AM
            [days, IA, IC]= unique(floor(start_times/1000/3600/24 + session_times),'rows');
            
            % remove today
            [~, idx] = intersect(days,round(now));
            if ~isempty(idx); IC = IC(IC~=IA(idx)); days = days(IA~=IA(idx));end
            
            % fetch lick times
            [ltimes, lprobe] =(fetchn(beh.Lick & key,'time','probe','ORDER BY time'));
            
            % loop through all days of the session
            for iday = unique(IC)'
                
                % intialize
                tuple = key;
                tuple.day = datestr(days(iday),'YYYY-mm-dd HH:MM:SS');
                
                % seperate target and distructor trials
                TARGET_idx = rprobe(IC==iday)==1;
                DISTR_idx = rprobe(IC==iday)==2;
                TARGET_sum = sum(TARGET_idx);
                DISTR_sum = sum(DISTR_idx);
                trial_start_times = start_times(IC==iday);
                trial_end_times = end_times(IC==iday);
                rew_probe = rprobe(IC==iday);
                
                % disregard day with very few trials
                if (TARGET_sum + DISTR_sum) < 20; continue; end
                
                % calculate start and stop times
                restr_idx = ltimes>min(trial_start_times) & ltimes<max(trial_end_times);
                day_probe = lprobe(restr_idx);
                day_lick_times = ltimes(restr_idx);
                day_resp = zeros(length(trial_start_times),1);
                for itrial = 1:length(trial_start_times)
                    %                    day_resp(itrial) = any(day_lick_times>=trial_start_times(itrial) & day_lick_times<trial_end_times(itrial));
                    resp = day_probe(find(day_lick_times>=trial_start_times(itrial) & day_lick_times<trial_end_times(itrial),1,'first'));
                    %resp = lprobe(find(day_lick_times>=trial_start_times(itrial) & day_lick_times<trial_start_times(itrial)+170,1,'first'));

                    if ~isempty(resp)
                        day_resp(itrial) = resp;
                    end
                end
                
                % compute correct and wrong lick probabilities
                HT_rate = sum(day_resp==1 & rew_probe==1)/TARGET_sum; %  p(lick|target)
                FA_rate = sum(day_resp==1 & rew_probe==2)/DISTR_sum; %  p(lick|distractor)
                tuple.ht = sum(day_resp==1 & rew_probe==1);
                tuple.fa = sum(day_resp==1 & rew_probe==2);
                tuple.ms = TARGET_sum - sum(day_resp==1 & rew_probe==1);
                tuple.cr = DISTR_sum - sum(day_resp==1 & rew_probe==2);
                
                % bootstrap licks
                rResp = cell2mat(arrayfun(@(x) day_resp(randperm(end))',1:1000,'uni',0)');
                r_HIT = sum(rResp(:,TARGET_idx & rew_probe==1),2)>0./TARGET_sum;
                r_FA = sum(rResp(:,DISTR_idx & rew_probe==2),2)>0./DISTR_sum;
                tuple.pval = 1 - mean(r_HIT-r_FA<HT_rate-FA_rate);
                
                % calculate metrics
                zH = norminv(min([max([HT_rate 0.01^40]) 0.99^40]),0,1);
                zF = norminv(min([max([FA_rate 0.01^40]) 0.99^40]),0,1);
                vH = HT_rate * (1-HT_rate) / (TARGET_sum*exp(-zH*zH/2)^2/(2*pi));
                vF = FA_rate * (1-FA_rate) / (DISTR_sum*exp(-zF*zF/2)^2/(2*pi));
                tuple.dprime = (zH - zF)/sqrt(2);
                tuple.perf = (tuple.ht + tuple.cr)/(TARGET_sum + DISTR_sum);
                tuple.dconf = 1.96*sqrt(vH+vF)/sqrt(2);
                
                % insert
                insert(self,tuple)
            end
        end
    end
end