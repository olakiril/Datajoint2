%{
# Performance with GonoGo tasks for each day
->beh.Session
day                   : timestamp                      # day starts at 12AM
---
HT                    : mediumint                      # hit trials
FA                    : mediumint                      # false alarm trials
MS                    : mediumint                      # miss trials
CR                    : mediumint                      # correct rejection trials
pval                  : float                          # bootstrap p value of HIT_RATE > FA_RATE
dprime                : float                          # Sensitivity d? = z(Hit) ? z(FA)
perf                  : float                          # performance (HITS + CR/ALL TRIALS)
dconf                 : float                          # 95% d' Confidence
%}

classdef LickPerformance < dj.Computed
    properties
        keySource = beh.Session & (beh.RewardCond & 'probe=0') & (beh.RewardCond & 'probe>0') & (beh.MovieClipCond & ...
                'movie_name="obj1v5" OR movie_name="obj2v5" OR movie_name="obj3v5" OR movie_name="obj3v5"')
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % fetch session times
            [start_times, session_times, end_times, rprobe] = ...
                fetchn((beh.Trial & key) * beh.Session * (beh.Condition & key) * (beh.RewardCond & key),...
                'start_time','session_tmst','end_time','probe','ORDER BY start_time');
            session_times = datenum(session_times,'YYYY-mm-dd HH:MM:SS'); 
            
            % start day at 12AM
            [days, ~, IC]= unique(floor(start_times/1000/3600/24 + session_times),'rows');
            
            % fetch lick times
            ltimes =(fetchn(beh.Lick & key,'time','ORDER BY time'));
           
            % loop through all days of the session
            for iday = unique(IC)'
                
                % intialize
                tuple = key;
                tuple.day = datestr(days(iday),'YYYY-mm-dd HH:MM:SS');
                
                % seperate target and distructor trials
                TARGET_idx = rprobe(IC==iday)>0;
                DISTR_idx = rprobe(IC==iday)==0;
                TARGET_sum = sum(TARGET_idx);
                DISTR_sum = sum(DISTR_idx);
                day_start = start_times(IC==iday);
                day_end = end_times(IC==iday);
                
                % calculate start and stop times
                day_lick_times = ltimes(ltimes>min(day_start) & ltimes<max(day_end));
                day_resp = false(1,length(day_start));
                for itrial = 1:length(day_start)
                    day_resp(itrial) = any(day_lick_times>=day_start(itrial) & day_lick_times<day_end(itrial));
                end
                
                % compute correct and wrong lick probabilities
                HT_rate = sum(day_resp(TARGET_idx))/TARGET_sum; %  p(lick|target)
                FA_rate = sum(day_resp(DISTR_idx))/DISTR_sum; %  p(lick|distractor)
                tuple.HT = sum(day_resp(TARGET_idx)); 
                tuple.FA = sum(day_resp(DISTR_idx)); 
                tuple.MS = TARGET_sum - sum(day_resp(TARGET_idx)); 
                tuple.CR = DISTR_sum - sum(day_resp(DISTR_idx));
                
                % bootstrap licks
                rResp = cell2mat(arrayfun(@(x) day_resp(randperm(end)),1:1000,'uni',0)');
                r_HIT = sum(rResp(:,TARGET_idx),2)./TARGET_sum;
                r_FA = sum(rResp(:,DISTR_idx),2)./DISTR_sum;
                tuple.pval = 1 - mean(r_HIT-r_FA<HT_rate-FA_rate);
                
                % calculate metrics
                zH = norminv(min([max([HT_rate 0.01^40]) 0.99^40]),0,1);
                zF = norminv(min([max([FA_rate 0.01^40]) 0.99^40]),0,1);
                vH = HT_rate * (1-HT_rate) / (TARGET_sum*exp(-zH*zH/2)^2/(2*pi));
                vF = FA_rate * (1-FA_rate) / (DISTR_sum*exp(-zF*zF/2)^2/(2*pi));
                tuple.dprime = (zH - zF)/sqrt(2);
                tuple.perf = (tuple.HT + tuple.CR)/(TARGET_sum + DISTR_sum);
                tuple.dconf = 1.96*sqrt(vH+vF)/sqrt(2);
                
                % insert 
                insert(self,tuple)
            end
        end
    end
end