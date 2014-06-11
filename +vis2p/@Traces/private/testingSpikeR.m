%% get quality of patched cells

if 1==1 % manualy inspected cell location of the patched cell
sI = {};

sI{end+1,1} = '2012-07-29';
sI{end,2} =10;

sI{end+1,1} = '2012-07-29';
sI{end,2} =8;

sI{end+1,1} = '2012-07-28';
sI{end,2} =14;

sI{end+1,1} = '2012-07-28';
sI{end,2} =13;

sI{end+1,1} = '2012-07-28';
sI{end,2} =9;

sI{end+1,1} = '2012-07-28';
sI{end,2} =8;

sI{end+1,1} = '2012-07-24';
sI{end,2} =15;

sI{end+1,1} = '2012-06-29';
sI{end,2} =19;

sI{end+1,1} = '2012-06-29';
sI{end,2} =18;

sI{end+1,1} = '2012-06-13';
sI{end,2} =28;

sI{end+1,1} = '2012-06-13';
sI{end,2} =27;

sI{end+1,1} = '2012-06-13';
sI{end,2} =25;

sI{end+1,1} = '2012-06-13';
sI{end,2} =23;

sI{end+1,1} = '2012-06-13';
sI{end,2} =13;

sI{end+1,1} = '2012-06-14';
sI{end,2} =6;

sI{end+1,1} = '2012-07-28';
sI{end,2} =17;

sI{end+1,1} = '2012-07-28';
sI{end,2} =4;

sI{end+1,1} = '2012-07-28';
sI{end,2} =3;

sI{end+1,1} = '2012-06-28';
sI{end,2} =16;

sI{end+1,1} = '2012-06-28';
sI{end,2} =5;

sI{end+1,1} = '2012-06-27';
sI{end,2} =5;

sI{end+1,1} = '2012-06-19';
sI{end,2} =6;

sI{end+1,1} = '2012-06-14';
sI{end,2} =8;

sI{end+1,1} = '2012-06-14';
sI{end,2} =6;

sI{end+1,1} = '2012-06-14';
sI{end,2} =4;

sI{end+1,1} = '2012-06-13';
sI{end,2} =21;

sI{end+1,1} = '2012-06-13';
sI{end,2} =17;

sI{end+1,1} = '2012-06-13';
sI{end,2} =13;

sI{end+1,1} = '2012-06-13';
sI{end,2} = 11;

end

bin = 100;
trace_opt = 27;

obj = EphysTraces.*Scans('scan_prog = "MpScan"');
keys = fetch(obj);
icell = 0;
skeys = []; gcorr = []; acorr = [];
keyind= true(length(keys),1);
scans = [];
for ikey = 1:length(keys)
    key = keys(ikey);
    k = [];
    cellnum = fetch1(Scans(key),'cell_patch');
    
    % see whether there is an aod scan with the same cell
    k.exp_date = key.exp_date;
    k.cell_patch = cellnum;
    k.scan_prog = 'AOD';
    if isempty(Scans(k));keyind(ikey) = false;continue;end
    
    % get mpscan data
    gfps = fetchn(Movies(key),'fps');
    gtraces  = getCaTraces(EphysTraces(key),'trace_opt',trace_opt);
    gspikes = fetchn(Traces(['masknum = 0 and trace_opt = ' num2str(trace_opt)],key),'trace');
    
    icell = icell+1;
    
dat{icell} =  k.exp_date;
scan(icell) = key.scan_idx;
    gcorr(icell) = corr(trresize(gtraces{1}(2:end),gfps,bin,'binsum'),...
        trresize(gspikes{1}(1:end-1),gfps,bin,'binsum')); %#ok<*SAGROW>
mf(icell) = mean(gspikes{1});
key.trace_opt = trace_opt;
qual(icell) = nanmean(fetchn(Traces(key),'quality'));
end

mean(gcorr)

%%   % get traces and frames


if 1==1 % manualy inspected cell location of the patched cell
sI = {};

sI{end+1,1} = '2012-07-29';
sI{end,2} =10;

sI{end+1,1} = '2012-07-29';
sI{end,2} =8;

sI{end+1,1} = '2012-07-28';
sI{end,2} =14;

sI{end+1,1} = '2012-07-28';
sI{end,2} =13;

sI{end+1,1} = '2012-07-28';
sI{end,2} =9;

sI{end+1,1} = '2012-07-28';
sI{end,2} =8;

sI{end+1,1} = '2012-07-24';
sI{end,2} =15;

sI{end+1,1} = '2012-06-29';
sI{end,2} =19;

sI{end+1,1} = '2012-06-29';
sI{end,2} =18;

sI{end+1,1} = '2012-06-13';
sI{end,2} =28;

sI{end+1,1} = '2012-06-13';
sI{end,2} =27;

sI{end+1,1} = '2012-06-13';
sI{end,2} =25;

sI{end+1,1} = '2012-06-13';
sI{end,2} =23;

sI{end+1,1} = '2012-06-13';
sI{end,2} =13;

sI{end+1,1} = '2012-06-14';
sI{end,2} =6;

sI{end+1,1} = '2012-07-28';
sI{end,2} =17;

sI{end+1,1} = '2012-07-28';
sI{end,2} =4;

sI{end+1,1} = '2012-07-28';
sI{end,2} =3;

sI{end+1,1} = '2012-06-28';
sI{end,2} =16;

sI{end+1,1} = '2012-06-28';
sI{end,2} =5;

sI{end+1,1} = '2012-06-27';
sI{end,2} =5;

sI{end+1,1} = '2012-06-19';
sI{end,2} =6;

sI{end+1,1} = '2012-06-14';
sI{end,2} =8;

sI{end+1,1} = '2012-06-14';
sI{end,2} =6;

sI{end+1,1} = '2012-06-14';
sI{end,2} =4;

sI{end+1,1} = '2012-06-13';
sI{end,2} =21;

sI{end+1,1} = '2012-06-13';
sI{end,2} =17;

sI{end+1,1} = '2012-06-13';
sI{end,2} =13;

sI{end+1,1} = '2012-06-13';
sI{end,2} = 11;

end

bin = 100;
trace_opt = 21;
% ikey = 35;
key = [];
key.exp_date = '2012-07-24';
key.scan_idx = 13;
% obj = EphysTraces.*Scans('scan_prog = "MpScan"');
% keys = fetch(obj);
icell = 0;
skeys = []; gcorr = []; acorr = [];
% keyind= true(length(keys),1);
scans = [];
% for ikey = 1:length(keys)
%     key = keys(ikey);

%     k = [];
%     cellnum = fetch1(Scans(key),'cell_patch');
%     
%     % see whether there is an aod scan with the same cell
%     k.exp_date = key.exp_date;
%     k.cell_patch = cellnum;
%     k.scan_prog = 'AOD';
%     if isempty(Scans(k));keyind(ikey) = false;continue;end
    
    clear traces
     [masknums,traces(:,:,1),traces(:,:,2),traces(:,:,3)] = ...
        fetchn(MaskTraces(key,'masknum > 0'),...
        'masknum','calcium_trace','red_trace','annulus_trace');
    fps = fetch1(Movies(key),'fps');
    opts = [3 17 22 21 24];
    TRACES = cell(length(opts),1);
    for iTopt = 1:length(opts)
        key.trace_opt = opts(iTopt);
        traceOpt = fetch(TracesOpt(key),'*'); 
        options = fetch1(TracesOpt(key),'trace_computation');
        options = strread(options, '%s','delimiter',',');
        
        % filter traces
        traceTypes = size(traces,3);
        tracesP = double([traces{:}]);  % traces are in columns
        tracesP = reshape(tracesP,size(tracesP,1),[],traceTypes);
        for iopt = 1:length(options)
            [tracesP, qual]= eval([options{iopt} '(tracesP,fps,traceOpt)']);
        end
        TRACES{iTopt} = tracesP;
    end
    
      
    % get mpscan data
    gfps = fetchn(Movies(key),'fps');
%     gtraces  = getCaTraces(EphysTraces(key),'trace_opt',trace_opt);
    key = rmfield(key,'trace_opt');
%     gspikes = fetchn(Traces(['masknum = 0 and trace_opt = 17'],key),'trace');
%     
    icell = icell+1;
    
     key.masknum =fetch1(EphysTraces(key),'masknum');
%     gcorr(icell) = corr(trresize(gtraces{1},gfps,bin,'binsum'),...
%         trresize(gspikes{1},gfps,bin,'binsum')); %#ok<*SAGROW>
% 
% end
% 
% mean(gcorr)

%%
close all
figure
s(1) = subplot(221);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{1},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{2},(1:size(TRACES{2},2))/2),'color',[1 0 0])
plot(normalize(gspikes{1})+key.masknum/2,'color','b')
title([num2str(opts(2)) ' ' num2str(corr(TRACES{2}(2:end,key.masknum),gspikes{1}(1:end-1)))])
set(gca,'ytick',[1:size(TRACES{1},2)]/2,'yticklabel',kurtosis(TRACES{1}));

s(2) = subplot(222);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{3},(1:size(TRACES{2},2))/2),'color',[1 0 0])
plot(normalize(gspikes{1})+key.masknum/2,'color','b')
title([num2str(opts(3)) ' ' num2str(corr(TRACES{3}(2:end,key.masknum),gspikes{1}(1:end-1)))])
set(gca,'ytick',[1:size(TRACES{1},2)]/2,'yticklabel',kurtosis(TRACES{1}));

s(3) = subplot(223);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{4},(1:size(TRACES{2},2))/2),'color',[1 0 0])
plot(normalize(gspikes{1})+key.masknum/2,'color','b')
title([num2str(opts(4)) ' ' num2str(corr(TRACES{4}(2:end,key.masknum),gspikes{1}(1:end-1)))])
set(gca,'ytick',[1:size(TRACES{1},2)]/2,'yticklabel',kurtosis(TRACES{3}));

s(4) = subplot(224);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{5},(1:size(TRACES{2},2))/2),'color',[1 0 0])
plot(normalize(gspikes{1})+key.masknum/2,'color','b')
title([num2str(opts(5)) ' ' ...
    num2str(corr(TRACES{5}(2:end,key.masknum),gspikes{1}(1:end-1)))])
set(gca,'ytick',[1:size(TRACES{1},2)]/2,'yticklabel',kurtosis(TRACES{3}));
title(num2str(opts(5)))
linkaxes(s)


   