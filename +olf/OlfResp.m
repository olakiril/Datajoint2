%{
olf.OlfResp (computed) #
-> olf.Traces
-> olf.RespOpt
---
response                    : mediumblob                    # response matrix [stimuli trials]
responsiveness              : double                        # ttest p value between all ON and OFF periods
stimuli                     : mediumblob                    # stimuli
%}


classdef OlfResp < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel = (olf.Traces*olf.RespOpt('process = "yes"'))
    end
    
    
    methods
        function self = OlfResp(varargin)
            self.restrict(varargin{:})
        end
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
            
            trace = fetch1(olf.Traces & key,'trace');
            [stimTrials,fps,nslices] = fetch1(olf.Scan & key,'trials','fps','zsize');
            stimTrials = stimTrials(1:nslices:end);
            [stims, trials] = fetchn(olf.StimPeriods & key,'stimulus','trial');
            [on,off] = fetchn(olf.RespOpt & key,'response_period','baseline_period');
            ustims = unique(stimTrials);
            mxtrial = max(ustims([1 diff(ustims)]==1));
            if mxtrial<0.8*length(stims)
                display('Too many trials missing!')
            end
            stims = stims(1:mxtrial);
            trials = trials(1:mxtrial);
            uniStims = unique(stims);
            
            R_ON = [];
            R_OFF = [];
            for iuni = 1:length(uniStims);
                stim = uniStims(iuni);
                uni_trials = trials(strcmp(stims,stim));
                for itrial = 1:length(uni_trials)
                    tstart = find(stimTrials == uni_trials(itrial),1,'first');
                    if off
                        OFF = mean(trace(tstart-round(fps*off/1000):tstart-1));
                    else
                        OFF = 0;
                    end
                    R_ON{iuni,itrial} = mean(trace(tstart:tstart+round(fps*on/1000)-1));
                    R_OFF{iuni,itrial} = OFF;
                end
            end
            
            index = ~any(cellfun(@isempty,R_ON));
            
            tuple = key;
            ron = cell2mat(R_ON(:,index));
            roff = cell2mat(R_OFF(:,index));
            
            [~,tuple.responsiveness] = ttest(ron(:),roff(:));
            tuple.response = ron - roff;
            tuple.stimuli = uniStims;
            
            insert( obj, tuple );
            
        end
    end
    
    
    methods
        function [resp,stim] = measureResponse(obj)
            keys = fetch(obj);
            for ikey = 1:length(keys)
                [R,uniStims] = fetch1(obj & keys(ikey),'response','stimuli');
                [st_idx,st_name] = fetchn(olf.Stimuli & keys(ikey),'stimulus_index','stimulus_description');
                mix = st_idx(~cellfun(@isempty,strfind(st_name,'Mix')));
                cof = st_idx(~cellfun(@isempty,strfind(st_name,'Coffee')));
                cof = cellfun(@str2num,uniStims)==cof;
                mix = cellfun(@str2num,uniStims)== mix;
                stim = st_name(1:size(R,1));
                t = ttest(R');
                if t(cof)||t(mix)
                    resp(:,ikey) = t;
                else
                    resp(:,ikey) = nan(size(t));
                end
                
            end
        end
        
        function saveData(obj)
            
            A = [];
            k = [];
            k.trace_opt = 3;
            
            A{1,1} = 'Mouse #';
            A{2,1} = 'Exp Date';
            A{3,1} = 'Neurons Recorded';
            A{4,1} = '% Responsive';
            A{5,1} = 'From the responsive:';
            A{6,1} = 'Pentanol';
            A{7,1} = 'Anisole';
            A{8,1} = 'Acetephenone';
            A{9,1} = 'Isoamyl Acetate';
            mice = unique(fetchn(olf.Mice & obj,'mouse_id'));
            for imouse = 1:length(mice);
                k.mouse_id = imouse;
                
                exp_date = fetchn(olf.Session & obj & k,'session_timestamp');
                [r,s] = measureResponse(olf.OlfResp & obj & k);
                R =  mean(~isnan(r'));
                S = nanmean(r');
                A{1,imouse+1} = mice(imouse);
                A{2,imouse+1} = exp_date{1}(1:10);
                A{3,imouse+1} = size(r,2);
                A{4,imouse+1} = R(1)*100;
                A{6,imouse+1} = S(1)*100;
                A{7,imouse+1} = S(2)*100;
                A{8,imouse+1} = S(3)*100;
                A{9,imouse+1} = S(4)*100;
                
            end
            
            
            filename = 'OlfResults.xlsx';
            sheet = 1;
            xlswrite(filename,A,sheet)
            
        end
        
        function saveRawData(obj)
                      
            Data = [];
            mice = unique(fetchn(olf.Mice & obj,'mouse_id'));
            for imouse = 1:length(mice);
                k = [];
                k.trace_opt = 3;
                k.mouse_id = imouse;
                
                [sessions,depth] = fetchn(olf.Session & obj & k,'session_timestamp','depth');

                Data.mice(imouse).mouse_id = imouse;
                Data.mice(imouse).exp_date = sessions{1}(1:10);
                
                for isession = 1:length(sessions)
                    key = k;
                    key.session_timestamp = sessions{isession};
                    Data.mice(imouse).sites(isession).time =  sessions{isession}(12:end);
                    Data.mice(imouse).sites(isession).depth = depth(isession);
                    
                    [R,uniStims,masknum] = fetchn(obj & key,'response','stimuli','masknum');
                    [sidx,sname] = fetchn(olf.Stimuli & key,'stimulus_index','stimulus_description');
                    stims = [];
                    for istim = 1:length(uniStims{1});
                       stims{istim} = sname{str2num(uniStims{1}{istim})==sidx};
                    end
                    for icell = 1:length(masknum)                      
                        Data.mice(imouse).sites(isession).cells(icell).masknum = masknum(icell);
                        Data.mice(imouse).sites(isession).cells(icell).stimuli = stims';
                        Data.mice(imouse).sites(isession).cells(icell).Response = R{icell};
                    end
                end
            end
            save('OlfRawData','Data')
        end
    end
end