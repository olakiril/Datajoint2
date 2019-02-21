%{
# miniature event detection
-> psp.Traces
-> psp.EventOpt
---
event_count           : int                       # performance
avg_amplitude         : float                     # average amplitude of psps
avg_rise_time         : float                     # average amplitude of psps
avg_decay_time        : float                     # average amplitude of psps
amplitude             : mediumblob                # amplitude of psp
rise_time             : mediumblob                # rise time
offset                : mediumblob                # x-offset of alpha function
decay_time            : mediumblob                # decay time
%}

classdef Events < dj.Computed
    
    properties
        keySource  = psp.EventOpt * psp.Traces
    end
    
    methods(Access=protected)
        function makeTuples(self, key)
            
            % fetch data
            [trace,trace_std] = fetch1(psp.Traces & key,'trace','std');
            fs = fetch1(psp.Exp & key,'fs');
            
            % fetch extraction parameters
            [lp,hp,amp_std, p.maxAmp, p.minTau, p.maxTau, p.minYOffset, p.maxYOffset, p.minDecay, p.maxDecay,p.derThresh, ...
                p.closestEPSPs, p.errThresh, p.dataFilterType, p.derFilterType, p.dataFilterLength, p.derFilterLength] = ...
                fetch1(psp.EventOpt & key,...
                'low_pass','high_pass','min_amp','max_amp','min_tau','max_tau','min_y_offset','max_y_offset','min_decay',...
                'max_decay','der_thresh','closestepsps','errthresh','datafiltertype','derfiltertype','datafilterlength','derfilterlength');
            p.minAmp = amp_std*trace_std;
            p.debugging=0;
            p.dataStart= 1;
            p.forceDisplay = 0;
            p.noFit=0;
            
            % band pass filter
            if trace_std>10^-13
                ftrace = bandpass(trace,[hp lp],fs,'ImpulseResponse','iir');
            else
                ftrace = trace;
            end
            
            % detect minis
            PSPs = psp.detectPSPs(-double(ftrace'),0,p);
            
            % insert
            key.amplitude = PSPs(:,1);
            key.rise_time = PSPs(:,2);
            key.offset = PSPs(:,3);
            key.decay_time = PSPs(:,4);
            key.avg_amplitude = mean(PSPs(:,1));
            key.avg_rise_time = mean(PSPs(:,2));
            key.avg_decay_time = mean(PSPs(:,4));
            if  size(PSPs,1)==1 &&  PSPs(1)==0
                key.event_count = 0;
            else
                key.event_count = size(PSPs,1);
            end
            self.insert(key)
        end
    end
    
    methods
        function plot(self,varargin)
            params.downsample = 10;

            params = getParams(params,varargin);
            
            keys = fetch(self);
            for key = keys'
                trace = fetch1(psp.Traces & key,'trace');
                fs = fetch1(psp.Exp & key,'fs');
                [lp,hp] = fetch1(psp.EventOpt & key,'low_pass','high_pass');

                [amplitude,offset] = fetch1(psp.Events & key,'amplitude','offset');

                 % band pass filter
                ftrace = bandpass(trace(1:params.downsample:end),[hp lp],fs/params.downsample,'ImpulseResponse','iir');

                % plot
                figure
                plot((1:length(ftrace))*params.downsample/fs,-ftrace)
                hold on
                idx = amplitude~=0;
                if sum(idx)>0
                    plot(offset(idx)/fs,amplitude(idx),'*r','markersize',10)
                end
                name = sprintf('%s Probe:%d Trial:%d',key.parent,key.probe,key.trial);
                set(gcf,'name',name)
            end
        end
        
        function plotProbe(self,varargin)
            params.downsample = 10;

            params = getParams(params,varargin);
            
            keys = fetch(psp.Probes & self);
            for tkey = keys'
                trace = fetchn(psp.Traces & tkey,'trace');
                fs = fetch1(psp.Exp & tkey,'fs');
                [lp,hp] = fetchn(psp.EventOpt & tkey,'low_pass','high_pass');

                [amplitude,offset] = fetchn(psp.Events & tkey,'amplitude','offset');
                Traces = cell(1,length(amplitude));
                for itrial = 1:length(amplitude)
                    % band pass filter
                    Traces{itrial} = bandpass(trace{itrial}(1:params.downsample:end),[hp lp],fs/params.downsample,'ImpulseResponse','iir');
                    if itrial>1
                        offset{itrial} = offset{itrial} + length(trace{itrial-1})*(itrial-1);
                    end
                end
                ftrace = cell2mat(Traces');
                amplitude = cell2mat(amplitude);
                offset = cell2mat(offset);
                
                % plot
                figure
                plot((1:length(ftrace))*params.downsample/fs,-ftrace)
                hold on
                idx = amplitude~=0;
                if sum(idx)>0
                    plot(offset(idx)/fs,amplitude(idx),'*r','markersize',10)
                end
                name = sprintf('%s Probe:%d',tkey.parent,tkey.probe);
                set(gcf,'name',name)
            end
        end
    end
end