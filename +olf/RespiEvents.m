%{
#
-> olf.RespiRaw
-> olf.RespiOpt
---
peaks               : mediumblob               # frames per second of the recording
troughs             : mediumblob               # frames per second of the recording
peak_times          : mediumblob               # total frames recorded
trough_times        : mediumblob               # final sampling rate of the respiration signal
filtered_trace      : mediumblob               # raw respiration signal
%}


classdef RespiEvents < dj.Imported
    
    properties
        keySource = olf.RespiRaw * olf.RespiOpt
    end
    
    methods (Access=protected)
        function makeTuples(self, key)
            
            % get Parameters
            [high_pass_fr,...
                low_pass_fr,...
                peak_min_interval,...
                trough_min_interval,...
                peak_thr_prctile,...
                trough_thr_prctile,...
                peak_thr_factor,...
                trough_thr_factor,...
                peak_height_thr_prctile,...
                trough_height_thr_prctile,...
                peak_height_thr_factor,...
                trough_height_thr_factor] = fetch1(olf.RespiOpt & key,...
                'high_pass_fr',...
                'low_pass_fr',...
                'peak_min_interval',...
                'trough_min_interval',...
                'peak_thr_prctile',...
                'trough_thr_prctile',...
                'peak_thr_factor',...
                'trough_thr_factor',...
                'peak_height_thr_prctile',...
                'trough_height_thr_prctile',...
                'peak_height_thr_factor',...
                'trough_height_thr_factor');
            
            % get data
            [fps, data] = fetch1(olf.RespiRaw & key,'trace_fs','trace');
            
            % filter raw trace
            hpd = highpass(data,high_pass_fr,fps);
            lpd = lowpass(hpd,low_pass_fr,fps);
            
            %             % detect peaks & troughs
            %             peaks = ne7.dsp.spaced_max(lpd, peak_min_interval*fps);
            %             thr =  peak_thr_factor*quantile(lpd(peaks),peak_thr_quantile);
            %             peaks = peaks(lpd(peaks) > thr);
            %             troughs = ne7.dsp.spaced_max(-lpd, trough_min_interval*fps);
            %             thr =  trough_thr_factor*quantile(-lpd(peaks),trough_thr_quantile);
            %             troughs = troughs(-lpd(troughs) >thr);
            %
            %             % if multiple troughs are found without an intermediate peak, select only the last before a peak
            %             v = bsxfun(@minus,troughs,peaks');
            %             v(v<0) = inf;
            %             mnv = min(v);
            %             d = diff(troughs);
            %             t = troughs(1:end-1);
            %             t(d<mnv(2:end)) = [];
            %
            %             % if multiple peaks are found without an intermediate trough, select only
            %             % the first after a trough
            %             v = bsxfun(@minus,peaks,t');
            %             v(v<0) = inf;
            %             mnv = min(v);
            %             d = diff(peaks);
            %             p = peaks(2:end);
            %             p(d<mnv(2:end)) = [];
            %
            base = min(lpd);
            data = data-base;
            lpd = lpd-base;
            peaks = ne7.dsp.spaced_max(lpd, peak_min_interval*fps);
            troughs = ne7.dsp.spaced_max(-lpd, trough_min_interval*fps);
            
            % select troughs with a given height difference from peaks or
            % below a selected threshold
            v = bsxfun(@minus,troughs,peaks');
            [~,mni] = min(abs(v));
            height_diff = lpd(peaks(mni))-lpd(troughs);
            thr1 = prctile(height_diff,trough_height_thr_prctile)*trough_height_thr_factor;
            thr2 = prctile(lpd(troughs),trough_thr_prctile)*trough_thr_factor;
            troughs = troughs(height_diff > thr1 | lpd(troughs)<thr2);
            
            % select peaks with a given height difference from troughs or
            % above a selected threshold
            v = bsxfun(@minus,troughs,peaks');
            [~,mni] = min(abs(v'));
            height_diff = lpd(peaks)-lpd(troughs(mni));
            thr1 = prctile(height_diff,peak_height_thr_prctile)*peak_height_thr_factor;
            thr2 = prctile(lpd(peaks),peak_thr_prctile)*peak_thr_factor;
            peaks = peaks(height_diff > thr1 | lpd(peaks)>thr2);
            
            % select min/max of vallays and mountains
            t_idx = false(length(troughs),length(troughs));
            itrough = 1;
            while itrough < length(troughs)
                p_idx = troughs<peaks(find(peaks>troughs(itrough),1,'first'));
                if isempty(p_idx);p_idx = true;end
                idx= troughs>troughs(itrough) & p_idx;
                if all(~idx);itrough = itrough+1;
                else
                    idx(itrough) = true;
                    t_idx(itrough,:) = lpd(troughs)' > min(lpd(troughs(idx))) & idx;
                    itrough = find(idx,1,'last')+1;
                end
            end
            t = troughs;
            t(any(t_idx)) = [];
            
            p_idx = false(length(peaks),length(peaks));
            ipeak = 1;
            while ipeak < length(peaks)
                t_idx = peaks<troughs(find(troughs>peaks(ipeak),1,'first'));
                if isempty(t_idx);t_idx = true;end
                idx= peaks>peaks(ipeak) & t_idx;
                if all(~idx); ipeak = ipeak+1;
                else
                    idx(ipeak) = true;
                    p_idx(ipeak,:) = lpd(peaks)' < max(lpd(peaks(idx))) & idx;
                    ipeak = find(idx,1,'last')+1;
                end
            end
            p = peaks;
            p(any(p_idx)) = [];
            
            % insert
            key.peaks = lpd(p);
            key.troughs = lpd(t);
            key.peak_times = p;
            key.trough_times = t;
            key.filtered_trace = lpd;
            insert(self, key);
        end
    end
    
    methods
        function plot(self)
            [fps, data] = fetch1(olf.RespiRaw & self,'trace_fs','trace');
            key = fetch(self);
            [lpd,p,t] = fetch1(self,'filtered_trace','peak_times','trough_times');
            trials = fetch1(olf.Sync & self,'trials');
            trial_fps = fetch1(reso.ScanInfo & key,'fps');
            data = data - mean(data);
            mn = mean(lpd);
            lpd = lpd - mn;
            
            % convert to voltage from int16
            lpd = lpd/(2^16/4)*2;
            data = data/(2^16/4)*2;
            
            %plot
            colors = cbrewer('qual','Set2',5);
            figure
            plot((1:length(data))/fps,data,'color',colors(3,:))
            hold on
            plot((1:length(lpd))/fps,lpd,'color',colors(2,:))
            plot(p/fps,ones(size(p))*prctile(lpd,70),'.','color',colors(1,:),'markersize',15)
            plot(t/fps,ones(size(t))*prctile(lpd,30),'.','color',colors(4,:),'markersize',15)
            yl = get(gca,'ylim');
            for itrial = 1:length(unique(trials))
                values = find(trials==itrial);
                plot(values/trial_fps,repmat(prctile(lpd,99),1,length(values)),'color',colors(5,:),'linewidth',3)
            end
            l = legend({'Raw data','Filtered Data','Expiration start','Inspiration start','stimulus'});
            set(l,'box','off')
            ylabel('Rel. temperature (Volts)')
            xlabel('Time(sec)')
            title('Breathing')
            set(gcf,'name',sprintf('Animal %d session %d scan %d breathing',key.animal_id,key.session,key.scan_idx))
            set(gca,'box','off')
        end
    end
end