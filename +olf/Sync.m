%{
olf.Sync (imported) #
-> olf.Session
animal_id         : int               # id (internal to database)
session           : smallint          # session index for the mouse
scan_idx          : smallint          # number of TIFF stack file
---
nframes           : int               # number of frames
ntrials           : int               # number oftrials
trials            : mediumblob        # decoded trial time
%}


classdef Sync < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = (olf.Session)
    end
    
    methods (Access=protected)
        function makeTuples(obj,key)
            
            [file, path, dur,it] = fetch1(olf.Session(key),'file_name','path','pulse_duration','intertrial_interval');
            mxtrial = max(fetchn(olf.StimPeriods & key,'trial'));
            filetype = getLocalPath(fullfile(path,sprintf('%s%s',file,'*.tif')));
            display(['Reading file: ' filetype])
            k = [];
            k.animal_id = key.mouse_id;
            k.filename = file;
            [key.session, key.scan_idx, key.animal_id] = ...
                fetch1(experiment.Scan & k,'session','scan_idx','animal_id');
            
            reader = ne7.scanimage.Reader5(filetype);
            sz = size(reader);
            
            % get Data
            data_phd = reader(:,:,sz(3),:,:);
            
            % filter photodiode
            key.trials = findTrials(obj,data_phd,reader.fps,dur,sz(4),mxtrial,it);
            key.trials = key.trials(1:sz(4):end);
            
            % fill values
            key.nframes = sz(5);
            key.ntrials = max(key.trials);
            
            % check consistency of trials
            nt = key.ntrials - (length(unique(key.trials))-1);
            if key.nframes<1000
                error(sprintf('Not enough frames'))
            elseif nt ~= 0
                error(sprintf('Missing %d trials',nt))
            end
            
            % insert Scan
            insert( obj, key );
            
        end
    end
    
    methods
        function self = Scan(varargin)
            self.restrict(varargin{:})
        end
        
        function trials = findTrials(obj,data_phd,fps,dur,nslices,mxtrial,it)
            
            phd = reshape(mean(data_phd,2),[],1);
            phd = (mean(phd)+min(phd))/2-phd > 0;
            if it==0
                flipdur = nslices*dur*fps*size(data_phd,1)/33/1000; % in line time
                prod = 2 .^ (0 : 15 )';
                stimtimes = 1:size(phd);
                stim = zeros(size(stimtimes));
                b = phd(round(flipdur*2/3+flipdur*(1:2:32) + (1:floor(length(phd)-flipdur*33))'));
                t = b*prod;
                tt = [];
                tt(1)= find(t==1,1,'first');
                stim(stimtimes > tt(1) & stimtimes < tt(1)+33*flipdur) = 1;
                for i = 2:mxtrial
                    pt = false(size(t));
                    pt(tt(end):end) = true;
                    f = find(t==i &  pt,1,'first');
                    if isempty(f); break; end
                    tt(i) = f;
                    stim(stimtimes > tt(i) & stimtimes < tt(i)+33*flipdur) = i;
                end
                trials = stim(1:size(data_phd,1):end);
            else
                flipdur = nslices*dur*fps*size(data_phd,1)/33/1000; % in line time
                all_times = find(diff(phd));
                start_times = find(diff(phd)>0);
                start_times = start_times((start_times+(flipdur*32+flipdur/2))<length(phd));
                prod = 2 .^ (0 : 15 )';
                istart = 1;
                lookuptimes = all_times(istart)+flipdur/2+flipdur*(1:2:32);
                b = phd(round(lookuptimes))';
                trial = b*prod;
                next = find(start_times>lookuptimes(end)+flipdur*1/2,1,'first'); % next = find(start_times>lookuptimes(end)+flipdur*3/2,1,'first'); 2017-06-12
                stimtimes = 1:size(phd);
                stim = zeros(size(stimtimes));
                stim(stimtimes > all_times(istart) & stimtimes < all_times(istart)+33*flipdur) = trial;
                previous_trial = trial;
                while start_times(next) < all_times(end)
                    nextlookuptimes = all_times(start_times(next)==all_times)+flipdur/2+flipdur*(1:2:32);
                    b = phd(round(nextlookuptimes))';
                    next_trial = b*prod;
                    if next_trial ~= previous_trial + 1
                        next = next+1;
                    else
                        stim(stimtimes > all_times(start_times(next)==all_times) ...
                            & stimtimes < all_times(start_times(next)==all_times)+33*flipdur) = next_trial;
                        next = find(start_times>nextlookuptimes(end)+flipdur,1,'first'); % next = find(start_times>nextlookuptimes(end)+flipdur*2,1,'first'); 2017-06-12
                        previous_trial = next_trial;
                    end
                    if next>length(start_times)
                        break
                    end
                end
                trials = stim(1:size(data_phd,1):end);
            end
            
            %%
            %                         flipdur = nslices*dur*fps*size(data_phd,1)/33/1000; % in line time
            %                         all_times = find(diff(phd));
            %                         start_times = find(diff(phd)>0);
            %                         start_times = start_times((start_times+(flipdur*32+flipdur/2))<length(phd));
            %                         prod = 2 .^ (0 : 15 )';
            %                         istart = 1;
            %                         while istart <= length(all_times)
            %                             lookuptimes = all_times(istart)+flipdur/2+flipdur*(1:2:32);
            %                             b = phd(round(lookuptimes))';
            %                             trial = b*prod;
            %                             next = find(start_times>lookuptimes(end)+flipdur*1/2,1,'first'); % next = find(start_times>lookuptimes(end)+flipdur*3/2,1,'first'); 2017-06-12
            %                             stimtimes = 1:size(phd);
            %                             stim = zeros(size(stimtimes));
            %                             stim(stimtimes > all_times(istart) & stimtimes < all_times(istart)+33*flipdur) = trial;
            %                             previous_trial = trial;
            %                             while start_times(next) < all_times(end)
            %
            %                                 nextlookuptimes = all_times(start_times(next)==all_times)+flipdur/2+flipdur*(1:2:32);
            %                                 b = phd(round(nextlookuptimes))';
            %                                 next_trial = b*prod;
            %                                 if next_trial ~= previous_trial + 1
            %                                     break
            %                                 else
            %                                     istart = find(start_times(next)==all_times);
            %                                 end
            %                                 stim(stimtimes > all_times(start_times(next)==all_times) ...
            %                                     & stimtimes < all_times(start_times(next)==all_times)+33*flipdur) = next_trial;
            %                                 next = find(start_times>nextlookuptimes(end)+flipdur,1,'first'); % next = find(start_times>nextlookuptimes(end)+flipdur*2,1,'first'); 2017-06-12
            %                                 if start_times(next) > all_times(end)
            %                                     istart = length(all_times)+1;
            %                                     break
            %                                 end
            %                                 previous_trial = next_trial;
            %                             end
            %                             istart = istart+ 1;
            %                         end
            %                         trials = stim(1:size(data_phd,1):end);
            %%
        end
        
    end
end