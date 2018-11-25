%{
-> reso.FluorescenceTrace
-> olf.RespOpt
---
resp_on                    : mediumblob                    # on response matrix [stimuli trials]
resp_off                   : mediumblob                    # off response matrix [stimuli trials]
stimuli                     : mediumblob                    # stimuli
%}


classdef Responses < dj.Computed
    
    properties
         keySource = (reso.FluorescenceTrace*olf.RespOpt('process = "yes"') & olf.Sync)
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
            
            % fetch stuff
            trace = double(fetch1(reso.FluorescenceTrace & key,'trace'));
            fps = fetch1(reso.ScanInfo & key,'fps');
            stimTrials = fetch1(olf.Sync & key,'trials');
            [trials, stims] = fetchn( olf.StimPeriods & (olf.Sync & key),'trial','stimulus');
            [on,off,base, on_delay, off_delay, respiration] = fetchn(olf.RespOpt & key,...
                'response_period','off_response_period','baseline_period','response_delay','off_response_delay','respiration');
            
            % process traces
            hp = 0.02; 
            trace = trace + abs(min(trace(:)))+eps;
            trace = trace./ne7.dsp.convmirr(trace,hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
            trace = trace - prctile(trace,10);
            
            % restrict to inhalation periods
            if any(strcmp(respiration,{'peaks','troughs'}))
                % upsample traces to get presice bining
                trace= interp1(trace,linspace(1,length(trace),length(trace)*10),'linear');
                stimTrials= interp1(stimTrials,linspace(1,length(stimTrials),length(stimTrials)*10),'linear');
                fps = fps*10;
                [peak_times,trough_times,trace_fs] = fetch1(olf.RespiEvents * olf.RespiRaw  & key & 'respi_opt = 1','peak_times','trough_times','trace_fs');
                mn = min([length(trough_times) length(peak_times)]);
                rTimes = {trough_times(1:mn)/trace_fs*fps,peak_times(1:mn)/trace_fs*fps}; 
                if strcmp(respiration,'troughs'); rTimes = fliplr(rTimes);end
                if min(rTimes{1})>min(rTimes{2})
                     X1 = [rTimes{1}(1:end-1);rTimes{2}(1:end-1)];
                     X2 = [rTimes{1}(2:end);rTimes{2}(1:end-1)];
                else
                     X1 = [rTimes{1}(1:end-1);rTimes{2}(2:end)];
                     X2 = [rTimes{1}(2:end);rTimes{2}(2:end)];  
                end
                X1 = max(min(round(mean(X1)),length(trace)),1);
                X2 = max(min(round(mean(X2)),length(trace)),1);
                for i = 1:length(X1)
                    trace(X1(i):X2(i)) = nan;
                end
            end
            
            % compute stimuli
            ustims = unique(stimTrials);
            mxtrial = max(ustims([1 diff(ustims)]==1));
            if mxtrial<0.8*length(stims)
                disp('Too many trials missing!')
            end
            stims = stims(1:mxtrial);
            trials = trials(1:mxtrial);
            uniStims = unique(stims);
            
            % calculate responses
            R_ON = []; 
            R_OFF = [];
            for iuni = 1:length(uniStims)
                stim = uniStims(iuni);
                uni_trials = trials(strcmp(stims,stim));
                for itrial = 1:length(uni_trials)
                    tstart = find(stimTrials == uni_trials(itrial),1,'first');
                    tend = find(stimTrials == uni_trials(itrial),1,'last')+1;
                    if tend+round(fps*(off+off_delay)/1000)-1 > length(trace)
                        break
                    end
                    if base
                        ON_base = nanmean(trace(max([tstart-round(fps*base/1000) 1]):tstart-1));
%                        OFF_base = mean(trace(max([tend-round(fps*base/1000) 1]):tend-1));
                         OFF_base = ON_base;
                    else
                        ON_base = 0 ;
                        OFF_base = 0 ;
                    end
                    R_ON{iuni,itrial} = nanmean(trace(tstart:tstart+round(fps*(on+on_delay)/1000)-1)) - ON_base;
                    R_OFF{iuni,itrial} = nanmean(trace(tend:tend+round(fps*(off+off_delay)/1000)-1)) - OFF_base;
                end
            end
            
            % remove incomplete trials
            index = ~any(cellfun(@isempty,R_ON));
            
            % insert
            tuple = key;
            tuple.resp_on = cell2mat(R_ON(:,index));
            tuple.resp_off = cell2mat(R_OFF(:,index));
            tuple.stimuli = uniStims;
            insert( obj, tuple );
            
        end
    end
    
    methods (Static)
        function trace = dfof(trace,fps)
            hp = 0.02; 
            trace = trace + abs(min(trace(:)))+eps;
            trace = trace./ne7.dsp.convmirr(trace,hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
            trace = trace - prctile(trace,10); 
        end
    end
end