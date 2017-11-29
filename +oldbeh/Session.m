%{
beh.Session (manual) #
-> vis2p.Mice
session_timestamp: timestamp             #
---
trial_interval=null         : smallint                      # time before the next trial (sec)
response_interval=null      : smallint                      # time of unique response (ms)
response_period=null        : smallint                      # time to wait for response (sec)
bad_delay=null              : smallint                      # time to punish (sec)
exp_type=null               : enum('Images','Orientation','Calibrate','Freerun','BW') #
stimuli=null                : char(255)                     #
rewarded_stimuli=null       : char(255)                     #
setup=null                  : tinyint                       # probe number
%}


classdef Session < dj.Relvar
    
    methods
        function self = Session(varargin)
            self.restrict(varargin{:})
        end
        
        function plotLicks(self)
            screensize = get( 0, 'Screensize' );
            
            for itype = 1:2
                
                figure
                set(gcf,'position',[(itype-1)*screensize(3)/2 30 screensize(3)/2 screensize(4)*0.9]);
                if itype == 1;
                    set(gcf,'name','Correct Stimuli')
                else
                    set(gcf,'name','Incorrect Stimuli')
                end
                k = [];
                
                mice = unique(fetchn(self,'mouse_id'));
                for imouse = 1:length(mice)
                    
                    k.mouse_id = mice(imouse);
                    
                    sessions = fetch(self & k & beh.Licks & beh.StimPeriods & 'exp_type > "Freerun"');
                    if isempty(sessions);continue;end
                    times = fetchn(beh.Session & sessions,'session_timestamp');
                    sessions = sessions([ diff(datenum(times,'yyyy-mm-dd HH:MM:SS'))*60*24>1;true]);
                    
                    for isession = 1:length(sessions)
                        session = beh.Session & sessions(isession);
                        p_types = unique(fetchn(session,'rewarded_stimuli'));
                        p_names = unique(fetchn(session & k ,'stimuli'));
                        resp_period = fetch1(session & k,'response_period');
                        
                        if isempty(strfind(p_types,'.png'))
                            A = textread(['Z:\users\Manolis\Labview\Water Conditioning\stimuli\' p_types{1} '.txt'],'%s');
                            rperiods = A{2};
                            p_names =A{1};
                            if itype == 1
                                if strcmp(rperiods(1),' ');rperiods = rperiods(2:end);end
                                periods =  strsplit(rperiods,',');
                            else
                                periods = setxor(strsplit(rperiods,','),strsplit(p_names,','));
                            end
                        else
                            if itype == 1
                                rperiods = [p_types{:}];
                                if strcmp(rperiods(1),' ');rperiods = rperiods(2:end);end
                                periods =  strsplit(rperiods,',');
                            else
                                periods = setxor(strsplit(p_types{:},','),strsplit(p_names{:},','));
                            end
                        end
                        if isempty(periods);continue;end
                        tex = [];
                        for i = 1:length(periods)
                            tex = [tex sprintf('period_type = "%s"',periods{i})];
                            if i~=length(periods); tex = [ tex ' or ']; end
                        end
                        wtimes = double(fetchn(beh.StimPeriods  & fetch(session) & tex,'timestamp'));
                        ltimes = double(fetchn(beh.Licks & fetch(session),'timestamp'));
                        stime = min([wtimes(:);ltimes(:)]);
                        
                        wtimes = (wtimes-stime)/1000/60;
                        ltimes = (ltimes-stime)/1000/60;
                        subplot(length(sessions),length(mice),imouse+length(mice)*(isession-1))
                        for i = 1:length(wtimes);
                            if i ~=length(wtimes)
                                licks = ltimes(ltimes>wtimes(i) & ltimes<wtimes(i+1));
                            else
                                licks = ltimes(ltimes>wtimes(i));
                            end
                            
                            plot(licks-wtimes(i),ones(size(licks))*i,'.k')
                            hold on
                        end
                        set(gca,'ydir','reverse','box','off')
                        xlim([-resp_period/60 4*resp_period/60])
                        ylim([0 length(wtimes)+1])
                        plot([0 0],[1 length(wtimes)],'b')
                        plot([resp_period resp_period]/60,[1 length(wtimes)],'--r')
                        if imouse==1
                            ylabel('Trial #')
                        end
                        
                        if isession == length(sessions); xlabel('Time (min)');end
                        
                        if isession == 1
                            setup = fetch1(session,'setup');
                            title([num2str(mice(imouse)) ' ' num2str(setup)])
                        end
                    end
                end
            end
        end
    end
end