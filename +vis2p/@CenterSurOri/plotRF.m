function plotRF(obj,k,varargin)

params.background = 'stim';

params = getParams(params,varargin);

import vis2p.*

% get stims files
[kk.x, kk.y, kk.exp_date] = fetch1(vis2p.Scans(k),'x','y','exp_date');

stims = fetchn(VisStims('exp_type = "MouseDotMapping"').*vis2p.Scans(kk,'problem_type = "none!"'),'stim_file');
stimE = fetch1(VisStims(k,'exp_type = "CenterSurround"'),'stim_file');
kkk.masknum = k.masknum;
if length(stims)>1
    [dS, stidx, scan] = fetchn(vis2p.RFFit(kkk,'rf_opt_num = 3').*vis2p.Scans(kk,'problem_type = "none!"'),'snr','stim_idx','scan_idx');
    [~,idx] = max(dS);
    k.stim_idx = stidx(idx);
    k.scan_idx = scan(idx);
    stim = stims{idx};
else
    stim = stims{1};
end

% get params
if strcmp(params.background,'stim')
    
    monSize = stim.params.constants.monitorCenter*2;
    xnum = stim.params.constants.dotNumX;
    ynum = stim.params.constants.dotNumY;
    intRad = stimE.params.constants.internalRadius;
    extRad = stimE.params.constants.externalRadius;
    loc =  stimE.params.constants.location;
    
    % create stimulus
    [x,y] = meshgrid((1:monSize(1))-monSize(1)/2,(1:monSize(2))-monSize(2)/2);
    border = (x-loc(1)).^2 + (y-loc(2)).^2 <= intRad.^2 | ...
        (x-loc(1)).^2 + (y-loc(2)).^2 >= extRad.^2 ;
    
    imagesc(border)
    colormap([.5 .5 .5;1 1 1]);
    hold on
    axis image
    axis off
    hold on
    
    % plot RFs
    keys = fetch(vis2p.RFFit(k));
    
    for i = 1:length(keys)
        key = keys(i);
        plot(vis2p.RFFit(key),'pos',[xnum ynum monSize(1) monSize(2)])
    end
    
else
    % plot RFs
    keys = fetch(vis2p.RFFit(k));
    
    for i = 1:length(keys)
        key = keys(i);
        
        map = fetch1(RFMap(key),'onpoff_rf');
        imagesc(map)
        colormap gray
        hold on
        plot(vis2p.RFFit(key))
        axis off
        axis image
        
    end
    
end