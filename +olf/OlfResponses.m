%{
olf.OlfResponses (computed) #
-> preprocess.SpikesRateTrace
-> olf.RespOpt
---
resp_on                    : mediumblob                    # on response matrix [stimuli trials]
resp_off                   : mediumblob                    # off response matrix [stimuli trials]
stimuli                     : mediumblob                    # stimuli
%}


classdef OlfResponses < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel = (preprocess.SpikesRateTrace*olf.RespOpt('process = "yes"') & olf.Sync)
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
            
            % fetch stuff
            trace = fetch1(preprocess.SpikesRateTrace & key,'rate_trace');
            fps = fetch1(preprocess.PrepareGalvo & key,'fps');
            stimTrials = fetch1(olf.Sync & key,'trials');
            [trials, stims] = fetchn( olf.StimPeriods & (olf.Sync & key),'trial','stimulus');
            [on,off,base, on_delay, off_delay] = fetchn(olf.RespOpt & key,...
                'response_period','off_response_period','baseline_period','response_delay','off_response_delay');
            
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
                        ON_base = mean(trace(max([tstart-round(fps*base/1000) 1]):tstart-1));
%                        OFF_base = mean(trace(max([tend-round(fps*base/1000) 1]):tend-1));
                         OFF_base = ON_base;
                    else
                        ON_base = 0 ;
                        OFF_base = 0 ;
                    end
                    R_ON{iuni,itrial} = mean(trace(tstart:tstart+round(fps*(on+on_delay)/1000)-1)) - ON_base;
                    R_OFF{iuni,itrial} = mean(trace(tend:tend+round(fps*(off+off_delay)/1000)-1)) - OFF_base;
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
end