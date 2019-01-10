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
            %             hp = 0.02;
            %             trace = trace + abs(min(trace(:)))+eps;
            %             trace = trace./ne7.dsp.convmirr(trace,hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
            %             trace = trace - prctile(trace,10);
            trace = self.dfof(trace,fps);
            
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
    
    methods
        function [Data,Stims,keys,Resp] = getData(self,varargin)
            
            params.pre = 3;
            params.post = 11;
            params.bin = 100;
            params.mxtrial = 10;
            params.hp = 0.02;
            params.respiration = 0;
            
            params = getParams(params,varargin);
            
            [traces,fps,stimTrials, keys] = fetchn(reso.FluorescenceTrace * proj(reso.ScanInfo,'fps') * proj(olf.Sync,'trials') & self, 'trace','fps','trials');
            
            R = cell(1,length(keys));
            
            if nargout>3
                params.respiration = true;
            end
            
            parfor ikey = 1:length(keys)
                
                % fetch stuff
                key= keys(ikey);
                trace = double(traces{ikey});
                trace = olf.Responses.dfof(trace,fps(ikey),params.hp);
                [trials, stims] = fetchn( olf.StimPeriods & (olf.Sync & key),'trial','stimulus');
                
                % compute stimuli
                ustims = unique(stimTrials{ikey});
                uniStims = unique(stims);
                
                % calculate responses
                RRR = cellfun(@(x) permute(x,[1 3 2]),num2cell(nan(length(uniStims),params.mxtrial,round((params.pre+params.post)*1000/params.bin)),3),'uni',0);
                
                % get respiration if asked
                if params.respiration
                   [Respi,r_fps] = fetch1(olf.RespiRaw*olf.RespiEvents & key,'filtered_trace','trace_fs') ;
                   Respi = normalize(Respi);
                   respi = RRR;
                end
                
                for iuni = 1:length(uniStims)
                    stim = uniStims(iuni);
                    uni_trials = trials(strcmp(stims,stim) & trials<=max(ustims));
                    for itrial = 1:length(uni_trials)
                        stim_start = find(stimTrials{ikey} == uni_trials(itrial),1,'first');
                        tstart = stim_start - round(fps(ikey)*params.pre) - 1;
                        tend = stim_start + round(fps(ikey)*params.post) - 1;
                        if tend+round(fps(ikey)*params.post)-1 > length(trace)
                            continue
                        end
                        RRR{iuni,itrial} = interp1(trace,linspace(tstart,tend,round((params.pre+params.post)*1000/params.bin)),'linear');
                        
                        if params.respiration
                            respi{iuni,itrial} = interp1(Respi,linspace(tstart/fps(ikey)*r_fps,tend/fps(ikey)*r_fps,round((params.pre+params.post)*1000/params.bin)),'linear');
                        end
                    end
                end
                R{ikey} = RRR;
                Stims{ikey} = uniStims;
                if params.respiration
                    Resp{ikey} = respi;
                end
            end
            
            Data = permute(cell2mat(cellfun(@(x) permute(cell2mat(cellfun(@(xx) permute(xx,[1 3 2]),x(:,1:params.mxtrial),'uni',0)),[1 4 2 3]),R,'uni',0)),[4 2 1 3]); % [time cells stim trials]
            if params.respiration
                Resp = permute(cell2mat(cellfun(@(x) permute(cell2mat(cellfun(@(xx) permute(xx,[1 3 2]),x(:,1:params.mxtrial),'uni',0)),[1 4 2 3]),Resp,'uni',0)),[4 2 1 3]); % [time cells stim trials]
            end

        end
        
        function plot(self,varargin)
            
            params.bin = 4000;
            params.hp = 0.01;
            params.scale = 1;
            params.alpha = 0.4;
            params.trace_color = [0.5 0.5 0.5];
            params.respiration = 1;
            params.rscale = 2;
            params.prctile = 50;
            
            params = getParams(params,varargin);
            
            keys = fetch(olf.Sync & self);
           
            for ikey = 1:length(keys)
                 clf
                k = keys(ikey);
                
                [traces, tkeys] = fetchn(reso.FluorescenceTrace & k & self,'trace');
                tr_idx = floor(linspace(1,size(traces,1),30));
                traces = traces(tr_idx);
                tkeys = tkeys(tr_idx);
                
                fps = fetch1(reso.ScanInfo & k,'fps');
                ntraces= double(cell2mat(traces'))*params.scale;
                ntraces = olf.Responses.dfof(ntraces,fps,params.hp);
                d = max(1,round(params.bin/1000*fps));
                ntraces = convn(ntraces,ones(d,1)/d,'same');
                
                
%                 if params.prctile
%                     [~,sort_idx] = sort(nanmean(ntraces),'A');
%                     ntraces = ntraces(:,sort_idx);
%                     tkeys = tkeys(sort_idx);
%                     ntraces = ntraces(:,1:15);
%                     tkeys = tkeys(1:15);
%                 end
                hold on
                
                [~,sort_idx] = sort([tkeys.mask_id],'ascend');
                ntraces = ntraces(:,sort_idx);
                t = fetch1(olf.Sync & k,'trials');
                s = fetchn(olf.StimPeriods &  (olf.Sync & k),'stimulus');
                c_idx = cellfun(@(x) str2double(x),s);
                un_stim = unique(c_idx);
                
                colors = cbrewer('qual','Set2',length(un_stim));
                
                for itrial = 1:max(t)
                    time_idx = t==itrial;
                    patch([[1 1]*find(time_idx,1,'first') [1 1]*find(time_idx,1,'last')]/fps,...
                        [0 1 1 0]*size(ntraces,2)+1,colors(un_stim==c_idx(itrial),:),...
                        'facealpha',params.alpha,'edgecolor','none')
                end
                for itrace = 1:size(ntraces,2)
                    plot((1:size(ntraces,1))/fps,ntraces(:,itrace)+itrace,'color',params.trace_color)
                end
                
                if params.respiration
                   [respi,r_fps] = fetch1(olf.RespiRaw*olf.RespiEvents & k,'filtered_trace','trace_fs') ;
                   plot((0:length(respi)-1)/r_fps,normalize(respi)*params.rscale-0.5,'-','color',params.trace_color);
                end
                set(gcf,'name',sprintf('Animal:%d Session:%d Scan:%d',k.animal_id,k.session,k.scan_idx))
                shg
                if ikey<length(keys)
                    pause
                end
            end
            
        end
    end
    
    methods (Static)
        function trace = dfof(trace,fps,hp)
            if nargin<3 || isempty(hp)
                hp = 0.02;
            end
            trace = trace + abs(min(trace(:)))+eps;
            trace = trace./ne7.dsp.convmirr(trace,hamming(round(fps/hp)*2+1)/sum(hamming(round(fps/hp)*2+1)))-1;  %  dF/F where F is low pass
        end
        
        function stim_groups = findSimilar(stims, stim_classes)
            
            stim = cellfun(@str2num, stims,'uni',0);
            stim_idx = cellfun(@(x) find(any(strcmp(repmat(x,size(sd,1),1),repmat(sd,1,size(x,2)))')),stim_classes,'uni',0);
            
            cmb = combnk(1:size(stim,1),2);
            
            % create stimulus length index to select only equal stim
            stim_length = cellfun(@length,stim);
            length_idx = stim_length(cmb(:,1)) == stim_length(cmb(:,2));
            
            %
            cis = false(size(stim,1),1);
            trans = false(size(stim,1),1);
            for iclass = find(~cellfun(@isempty,stim_idx))
                is_class = cellfun(@(x) all(any(bsxfun(@rdivide,x,stim_idx{iclass}')==1)),stim,'uni',1);
                has_class = cellfun(@(x) any(any(bsxfun(@rdivide,x,stim_idx{iclass}')==1)),stim,'uni',1) & ~is_class;
                
                dif_idx = any(is_class(cmb)')' & length_idx & ~all(is_class(cmb)')' ;
                sam_idx =  (all(is_class(cmb)')' | all(has_class(cmb)')') & length_idx;
                
                cis(sam_idx) = true;
                trans(dif_idx) = true;
            end
        end
        
    end
end