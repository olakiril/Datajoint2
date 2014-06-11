function eye_dir = findEye(obj,key,eye_path)

% function eye_dir = findEye(obj,key,eye_path)
%
% Checks if any eye file overlaps and gets the parent directory of the file

% select eye global directory

import vis2p.*


if nargin<3; eye_path = 'M:/EyeCamera/';end
[frames, fps] = fetch1(Movies(key),'nframes','fps');
date_v = fetchTime(Scans(key));

% calculate scan period
ts = round(date_v(4)*3600 + date_v(5)*60 + date_v(6)); %timestamp in seconds
key.scan_idx = key.scan_idx+1;
if ~isempty(Scans(key));
    ndate_v = fetchTime(Scans(key));
    nts = round(ndate_v(4)*3600 + ndate_v(5)*60 + ndate_v(6)); %timestamp in seconds
    scanperiod = nts-ts;
else
    scanperiod = 3600*2; % seconds of waiting till I pressed the record button
end

% day's worth of seconds & mark scan time slot with ones
daysec = zeros(24*3600,1); 
daysec(ts:ts+scanperiod) = 1;

% find overlaping times
eye_dirs = dir([eye_path datestr(date_v,'yy-mm-dd') '*']);
sec_over = zeros(length(eye_dirs),1);
for idir = 1:length(eye_dirs)
    testsec = daysec;
    eye_v = datevec(eye_dirs(idir).name,'yy-mm-dd_HH-MM-SS');
    avi_inf = mmfileinfo([eye_path eye_dirs(idir).name '/eyemovie.avi']);
    avi_dur = round(avi_inf.Duration); % in seconds
    eye_ts = round(eye_v(4)*3600+eye_v(5)*60+eye_v(6));
    testsec(eye_ts:eye_ts+avi_dur) = testsec(eye_ts:eye_ts+avi_dur)+1;
    sec_over(idir) =  sum((testsec==2));
end
[over_m, over_i] = max(sec_over);

% check that there is at least 10% overlap
if over_m>frames*0.1/fps;
    eye_dir = [eye_path eye_dirs(over_i).name '/'];
else
    display 'No file found'
    eye_dir = [];
end

%%
ad  = zeros(length(eye_dirs),1);
for idir = 1:length(eye_dirs)
    testsec = daysec;
    eye_v = datevec(eye_dirs(idir).name,'yy-mm-dd_HH-MM-SS');
    avi_inf = mmfileinfo([eye_path eye_dirs(idir).name '/eyemovie.avi']);
    avi_dur = round(avi_inf.Duration); % in seconds
 ad(idir) = avi_dur;
end