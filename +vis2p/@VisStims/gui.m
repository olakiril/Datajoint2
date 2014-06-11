function gui( visStims )
assert( ~isempty( visStims ), 'Nothing to show.' );

%global variables
traceInfo.description = 'data for currently selected traces';
twoPhotonRGB = [];
cellList = [];
traceOptions = struct([]);
scanKeys = enumerate( visStims );
clear visStims;  % don't need it anymore

% GUI LAYOUT
mainFigure = figure(...
    'PaperUnits', 'inches',...
    'PaperPosition', [0 0 20 10],...
    'PaperSize', [20 10],...
    'Units','pixels',...
    'MenuBar','none',...
    'Toolbar','none',...
    'Position',[10 10 1350 780],...
    'Visible','off',...
    'Name', 'Scans gui',...
    'NumberTitle', 'off',...
    'Resize', 'off', 'WindowButtonUpFcn', {@pressCell} );

titleText = uicontrol(mainFigure,'Style','text','String','<no data loaded>','Units', 'pixels','Position',[20 750 900 20], 'FontSize', 14, 'fontweight', 'bold');
cellSelectionText  = uicontrol(mainFigure,'Style','text','String','<no cells selected>','Units', 'pixels','Position',[20 30 500 16], 'FontSize', 10, 'HorizontalAlignment', 'left', 'Callback', {@clearCellSelection});
plotCellsButton    = uicontrol(mainFigure,'Style','pushbutton','String','Plot selected cells',  'Position', [530 40 120 18], 'Callback', {@plotSelectedCells}, 'Enable', 'off' );
clearCellsButton   = uicontrol(mainFigure,'Style','pushbutton','String','Clear selected cells', 'Position', [530 20 120 18], 'Callback', {@clearCellSelection}, 'Enable', 'off' );

intrinsicAxes = axes( 'Parent', mainFigure, 'units', 'pixels', 'Position', [20, 730-530, 400, 400]);
twoPhotonAxes = axes( 'Parent', mainFigure, 'units', 'pixels', 'Position', [450, 720-512, 512, 512] );
histAxes      = axes( 'Parent', mainFigure, 'units', 'pixels', 'Position', [1000, 540, 300, 180]);
scatterAxes   = axes( 'Parent', mainFigure, 'units', 'pixels', 'Position', [1000, 230, 300, 200]);

editNote = uicontrol(mainFigure,'Style','edit',...
    'String', '',...
    'ToolTipString','Enter your comment here',...
    'Max',16,'Min',0,...
    'HorizontalAlignment', 'left',...
    'Position',[1025 105 300 60]);
editText           = uicontrol(mainFigure,'Style','text','String','Edit comment','Position',[1075 170 200,12]);
saveCommentButton  = uicontrol(mainFigure,'Style','pushbutton','String','Save comment', 'Position', [1200 80 110 18], 'Callback', {@saveComment}, 'Enable', 'on' );
printPdfButton     = uicontrol(mainFigure,'Style','pushbutton','String','Save PDF in curr dir', 'Position', [1200 50 120 18], 'Callback', {@printPdf}, 'Enable', 'on' );

scanStrings = {};
for key = scanKeys
    scanStrings{end+1} = sprintf('%s::%03d %s', key.exp_date, key.scan_idx, take(Experiments(key),'operator') );
end
key = scanKeys(1);

scanChoice = uicontrol(mainFigure ,'Style','listbox',...
    'String', scanStrings,...
    'Value',1,'Position',[10 620 180 110], 'Callback', {@goToNewScan} );

paramChoice = uicontrol(mainFigure ,'Style','popupmenu',...
    'String', '<empty>',...
    'Value',1,'Position',[720 45 350 20], 'Callback', {@refreshTwoPhoton} );

summaryText = uicontrol(mainFigure,'Style','text','String','<no data loaded>','Units', 'pixels','Position',[1000 445 260 64]...
    ,'FontSize', 12, 'HorizontalAlignment', 'left', 'FontWeight', 'demi');


% histogram menu
histChoices(1) = struct('code', 'pref oris' , 'display', 'Preferred Orientations (p<0.05)');
histChoices(2) = struct('code', 'pref dirs' , 'display', 'Preferred directions (p<0.05');
histChoices(3) = struct('code', 'skew'      , 'display', 'trace skewness' );
histChoices(4) = struct('code', 'vM SNR'    , 'display', 'von Mises SNRs' );
histChoices(5) = struct('code', 'eig evoked', 'display', 'evoked eigenvalues (takes 5 s)' );
histChoices(6) = struct('code', 'eig spont' , 'display', 'spont eigenvalues (takes 5 s)' );

histChoiceMenu = uicontrol( mainFigure, 'Style', 'popupmenu', 'String', {histChoices.display}, 'Units', 'pixels', 'Position', [1050 740 190 18], 'Callback', {@refreshHisto} );

% optical
opticalChoicePanel= uibuttongroup('Parent',mainFigure,'Title','Optical imaging choices','Units', 'pixels', 'Position',[100, 80, 195, 90],'SelectionChangeFcn', {@refreshIntrinsic} );
fluoroInsert = uicontrol(opticalChoicePanel,'Style','radiobutton','String','structural + fluo insert','Units','normalized','Position',    [.05 .8 .9 .2] );
opticChoice = uicontrol(opticalChoicePanel,'Style','radiobutton','String','structural map','Units','normalized','Position',    [.05 .6 .9 .2] );
elevChoice  = uicontrol(opticalChoicePanel,'Style','radiobutton','String','elevation map','Units','normalized','Position',    [.05 .4 .9 .2] );
azimChoice  = uicontrol(opticalChoicePanel,'Style','radiobutton','String','azimuth map','Units','normalized','Position',    [.05 .2 .9 .2] );

mapChoices = struct('code', {}, 'display', {} );
mapChoiceMenu = uicontrol( mainFigure, 'Style', 'popupmenu', 'Units', 'pixels', 'Position',[500, 720, 230, 20], 'Callback', {@refreshTwoPhoton} );
mapChoiceLabel = uicontrol( mainFigure, 'Style', 'text', 'String', 'Maps:', 'Units','pixels','Position',[450 720, 46,16]);

cellChoicePanel = uibuttongroup('Parent',mainFigure,'Title','Cell display','Units', 'pixels', 'Position',[700, 90, 295, 100]);
circCellChoice  = uicontrol(cellChoicePanel,'Style','checkbox','String','circle all detected cells (green dashes)'        ,'Units','normalized','value',0,'Position',[.05 .78 .9 .2], 'Callback', { @refreshTwoPhoton});
circRespChoice  = uicontrol(cellChoicePanel,'Style','checkbox','String','circle visually responsive cells (orange dashes)','Units','normalized','value',0,'Position',[.05 .59 .9 .2], 'Callback', { @refreshTwoPhoton});
barCellChoice      = uicontrol(cellChoicePanel,'Style','checkbox','String','orientation bar (p-value<0.05)'               ,'Units','normalized','value',0,'Position',[.05 .41 .9 .2], 'Callback', { @refreshTwoPhoton});
arrowCellChoice    = uicontrol(cellChoicePanel,'Style','checkbox','String','direction arrow (p-value<0.05)'                ,'Units','normalized','value',0,'Position',[.05 .23 .9 .2], 'Callback', { @refreshTwoPhoton});
hueCellChoice      = uicontrol(cellChoicePanel,'Style','checkbox','String','orientation hue (p-value<0.05)'               ,'Units','normalized','value',0,'Position',[.05 .05 .9 .2], 'Callback', { @refreshTwoPhoton});

goToNewScan

set( mainFigure, 'Visible', 'on' );




%% refreshIntrinsic callback
    function refreshIntrinsic( hObj, event )
        set( mainFigure, 'Pointer', 'watch' );

        switch 1
            case get( opticChoice, 'value' )
                img = take(Intrinsics(key),'structure_image');
                img = img';
                if ~isempty(img)
                    assert( all( size(img) == 512 ), 'Optical images should be 512x512' );
                    img = img - prctile( img(:), 1 );
                    img = min( 1, img/prctile( img(:), 30 )*0.4);
                    img = cat(3, img, img, img );
                end
            case get( fluoroInsert, 'value' )
                [img,cx,cy,scale,camAngle] = take( Intrinsics(key),'structure_image','cam_x','cam_y','cam_scale','cam_angle');
                img = img';
                if ~isempty(img) && ~isnan(cx)
                    sz = size( img );
                    assert( all(sz==512), 'Optical images should be 512x512'  );
                    img = img - prctile( img(:), 1 );
                    img = min( 1, img/prctile( img(:), 30 )*0.4);
                    [xi,yi] = meshgrid( ((0:sz(2)-1)/(sz(2)-1)-0.5)*scale, (0.5-(0:sz(1)-1)/(sz(1)-1))*scale );  %coordinates in um

                    [scanx,scany,scanRotation,greenImg,redImg,lens,mag] = take( Scans(key)...
                        ,'scanx','scany','scan_rotation','green_mean_img','red_mean_img','lens','mag');
                    xp = cos(camAngle)*xi - sin(camAngle)*yi + cx - scanx;
                    yp = sin(camAngle)*xi + cos(camAngle)*yi + cy - scany;
                    xr = cos( scanRotation )*xp - sin( scanRotation )*yp;
                    yr = sin( scanRotation )*xp + cos( scanRotation )*yp;
                    sz = size(greenImg);
                    [xi,yi] = meshgrid( ((0:sz(2)-1) - (sz(2)-1)/2)...
                        ,-((0:sz(1)-1) - (sz(1)-1)/2) );
                    xi = xi*15000/512/lens/mag;
                    yi = yi*15000/512/lens/mag;

                    greenImg = log(greenImg);
                    greenImg = greenImg-min(greenImg(:));
                    greenImg = greenImg./max(greenImg(:));
                    redImg = log(redImg);
                    redImg = redImg-min(redImg(:));
                    redImg = redImg./max(redImg(:));

                    zg = interp2( xi, yi, max(0,greenImg)*0.9+0.20, xr, yr, 'linear', 0 );
                    zr = interp2( xi, yi, max(0,redImg  )*0.9+0.20, xr, yr, 'linear', 0 );
                    img = cat( 3, img.*(zr==0) + zr, 0.9*img.*(zg == 0) + zg, img.*(zg==0) );
                end

            case get( elevChoice, 'value' )
                img = take( Intrinsics(key),'elevation_image')';
                if ~isempty( img )
                    assert( all(size(img)==512), 'Optical images should be 512x512' );
                    h = mod(angle(img)/2/pi+pi+0.32,1);
                    s = abs(img);
                    s = max(0, s/quantile( s(:),0.98));
                    img = hsv2rgb( cat(3, h, cat(3, s.^1.5, min(1,s ) ) ) );
                end
            case get( azimChoice, 'value' )
                img = take( Intrinsics(key), 'azimuth_image')';
                if ~isempty( img )
                    assert( all(size(img)==512), 'Optical images should be 512x512' );
                    h = mod(angle(img)/2/pi+pi+0.32,1);
                    s = abs(img);
                    s = max(0, s/quantile( s(:),0.98));
                    img = hsv2rgb( cat(3, h, cat(3, s.^1.5, min(1, s ) ) ) );
                end
        end
        if isempty( img )
            cla( intrinsicAxes );
        else
            axes( intrinsicAxes );
            imshow( img );
            xlim( [106 505] );
            ylim( [106 505] );
            %% show rectangle in intrinsic image
            if ~get( fluoroInsert, 'value' )
                [greenImg,scanRotation,x0,y0,lens,mag] = take( Scans(key)...
                    ,'green_mean_img','scan_rotation','scanx','scany','lens','mag');
                sz = size(greenImg);
                xp = [-1 -1  1  1 -1]/2*sz(2);
                yp = [-1  1  1 -1 -1]/2*sz(1);
                scanRotation =  scanRotation*pi/180;
                xi = cos( scanRotation )*xp + sin( scanRotation )*yp;
                yi =-sin( scanRotation )*xp + cos( scanRotation )*yp;

                xi = xi*15000/512/lens/mag + x0;  % cellx in manipulator coordinates
                yi = yi*15000/512/lens/mag + y0;  % celly in manipulator coordinates
                [camx,camy,camScale,camAngle] = take( Intrinsics(key)...
                    ,'cam_x','cam_y','cam_scale','cam_angle');
                xi = xi - camx;
                yi = yi - camy;
                xx = cos(camAngle)*xi + sin(camAngle)*yi;
                yy =-sin(camAngle)*xi + cos(camAngle)*yi;
                xx = xx/camScale;
                yy = yy/camScale;
                xx = xx*512 + 255.5;
                yy =-yy*512 + 255.5;
                hold on; plot(xx,yy,'Color', [1 1 0.3], 'LineWidth', 1); hold off;
            end
        end
        set( mainFigure, 'Pointer', 'arrow' );
    end




    function refreshHisto( hObj, event )

        cla( histAxes );
        axes( histAxes );
        switch histChoices( get( histChoiceMenu, 'Value' ) ).code;
            case 'pref oris'
                if ~isempty( traceInfo.oriSelective )
                    colors = hsv( 180 );
                    [counts, bins] = hist( traceInfo.pref_ori( traceInfo.oriSelective ), 7.5:15:180 );
                    for iBin = 1:length(counts)
                        hbar = bar( bins, counts.*(1:length(counts) == iBin), 1 );
                        set(hbar, 'FaceColor', colors(round(bins(iBin)),:) );
                        hold on;
                    end
                    hold off;
                    set( histAxes, 'XLim', [0 180], 'YLim', [0 max(counts)*1.1+1]);
                    xlabel( 'degrees orientation' );
                end

            case 'pref dirs'
                if ~isempty( traceInfo.dirSelective )
                    colors = hsv( 180 );
                    [counts,bins] = hist( traceInfo.pref_dir( traceInfo.dirSelective ), 7.5:15:360 );
                    for iBin = 1:length(counts)
                        hbar = bar( bins, counts.*(1:length(counts) == iBin), 1 );
                        set(hbar, 'FaceColor', colors(max(1,round(mod( bins(iBin), 180))),:) );
                        hold on;
                    end
                    hold off;
                    set( histAxes, 'XLim', [0 360], 'YLim', [0 max(counts)*1.1+1]);
                    xlabel( 'degrees direction' );
                end

            case 'skew'
                sk = take( CircStats(key), 'skewness' );
                if ~isempty( sk )
                    hist( sk, -1.4:0.1:3.0 );
                    xlim( [-1.44 3.06] );
                    xlabel( 'trace skewness' );
                end

            case 'vM SNR'
                vmSNR = take( CircStats(key), 'vm_snr');
                if ~isempty( vmSNR )
                    x=-2:0.1:0;
                    h = hist( log10(vmSNR), x );
                    bar( x, h );
                    xlim([-2.2 0.2]);
                    set( gca, 'XTick', x(1:10:end), 'XTickLabel', 10.^x(1:10:end) );
                    xlabel( 'von Mises SNR' );
                end

            otherwise
                error( 'Unknown histogram code');
        end
    end



%% refreshTwoPhoton callback
    function refreshTwoPhoton( hObj, event )
        set(mainFigure,'Pointer','watch');
        if ~isempty(traceOptions)
            key.trace_opt = traceOptions(get(paramChoice,'Value')).trace_opt;
        end
        [directions,trials] = take(StimPresents(key), 'direction','trial');
        [lens,mag,depth] = take(Scans(key),'lens','mag','depth');
        [pdays,anesthesia] = take(Experiments(key),'pdays','anesthesia');
        str=sprintf( 'Age p%3d. Anesth %s. Lens=%d. Mag=%1.2f. Depth = %d. %d directions. %d trials.'...
            ,pdays,anesthesia,lens,mag,depth,length(unique(directions)),length(unique(trials)));
        set(titleText,'String',str);
        twoPhotonRGB = [];

        %% make two-photon map
        code = mapChoices( get( mapChoiceMenu, 'Value' ) ).code;
        switch code
            case 'gray'
                img = take(Movies(key),'affine_green');
                if isempty(img)
                    img = take(Movies(key),'raw_green');
                end
                img = img - prctile( img(:), 0.2);
                img = img / prctile( img(:), 99.8);
                twoPhotonRGB = cat( 3, img, img, img );
            case 'red/green'
                [img,g] = take(Movies(key),'affine_red','affine_green');
                if isempty(img)
                    [img,g] = take(Movies(key),'raw_red','raw_green');
                end
                img = img - prctile( img(:), 0.2);
                img = img / prctile( img(:), 99.8)*0.6;
                g = g - prctile( g(:), 0.2);
                g = g / prctile( g(:), 99.8);
                twoPhotonRGB = cat( 3, img, g, zeros(size(g)) );

            otherwise
                if ~strncmp(code, 'cat',3)
                    error('Unknown mapMenu code');
                else
                    tmp = sscanf( code, 'cat %d ncomps %d' );
                    k = key;
                    k.cat = tmp(1);
                    k.ncomps = tmp(2);
                    [s,h] = take( TuningMaps(k), 'cos_pval_map','cos_phase_map');
                    h = h/180;
                    min(s(:));
                    s = s < 0.05;

                    v = take(Scans(key),'green_mean_img');
                    v = v - prctile( v(:), 0.2);
                    v = v / prctile( v(:), 99.8);
                    twoPhotonRGB = hsv2rgb( cat( 3, h, s, v) );
                end
        end

        cla( twoPhotonAxes );
        if ~isempty( twoPhotonRGB )
            axes( twoPhotonAxes );
            imshow( twoPhotonRGB );
        end

        % load trace information
        traceInfo = take(CircStats(key,'masknum>1')*CircCells,'masknum','img_x','img_y','cell_radius','vis_p','pref_ori','ori_p','pref_dir','dir_p');
        traceInfo.dirSelective = find([traceInfo.dir_p]<0.05);
        traceInfo.oriSelective = find([traceInfo.ori_p]<0.05);
        refreshHisto;


        %% update summary text
        nCells = length( traceInfo );
        set( summaryText, 'String',...
            sprintf( 'Detected %d cells.  \n%d orientationally selective (%3.2f%%).   \n%d directionally selective (%3.2f%%).' ...
            , nCells, length(traceInfo.oriSelective), 100*length(traceInfo.oriSelective)/nCells, length(traceInfo.dirSelective), 100*length(traceInfo.dirSelective)/nCells) );

        %% show cells
        axes( twoPhotonAxes );
        if get( hueCellChoice, 'value' ) && ~isempty( traceInfo.oriSelective )
            hold on;
            for iCell = traceInfo.oriSelective'
                rectangle( 'Position',[[traceInfo.img_x(iCell), traceInfo.img_y(iCell)] - traceInfo.cell_radius(iCell), traceInfo.cell_radius(iCell)*[2 2]]...
                    , 'EdgeColor', 'none', 'Curvature', [1 1], 'FaceColor', hsv2rgb( [traceInfo.pref_ori(iCell)/180, 1, 1] )  );
            end
            hold off;
        end

        if get( circCellChoice, 'value' )
            hold on;
            for iCell = 1:nCells
                rectangle( 'Position', [[traceInfo.img_x(iCell)  traceInfo.img_y(iCell)]-1.2*traceInfo.cell_radius(iCell), 1.2*traceInfo.cell_radius(iCell)*[2 2]], 'Curvature', [1,1]...
                    , 'EdgeColor', [0.0 0.7 0.0], 'LineWidth', 1);
            end
            hold off;
        end

        if get( circRespChoice, 'value' )
            hold on;
            for iCell = find(traceInfo.vis_resp_pvalue' < 0.05)
                hrect = rectangle( 'Position', [[traceInfo.img_x(iCell)  traceInfo.img_y(iCell)]-1.2*traceInfo.cell_radius(iCell), 1.2*traceInfo.cell_radius(iCell)*[2 2]], 'Curvature', [1,1]...
                    , 'EdgeColor', [1.0 0.7 0.0], 'LineWidth', 1);
            end
            hold off;
        end

        if get( barCellChoice, 'value' )
            hold on;
            for iCell = traceInfo.oriSelective'
                plot(traceInfo.img_x(iCell) + cos( traceInfo.pref_ori(iCell)*pi/180 )*[-1 1]*1.0*traceInfo.cell_radius(iCell),...
                     traceInfo.img_y(iCell) + sin( traceInfo.pref_ori(iCell)*pi/180 )*[-1 1]*1.0*traceInfo.cell_radius(iCell), 'Color', [0 0.1 0], 'LineWidth', 1 );
            end
            hold off;
        end

        if get( arrowCellChoice, 'value' )
            hold on;
            for iCell = traceInfo.dirSelective'
                plot( traceInfo.img_x(iCell), traceInfo.img_y(iCell), 'k.');
                plot( traceInfo.img_x(iCell) + [0, -sin( traceInfo.pref_ori(iCell)*pi/180 )*traceInfo.cell_radius(iCell)],...
                    traceInfo.img_y(iCell) + [0, -cos( traceInfo.pref_ori(iCell)*pi/180 )*traceInfo.cell_radius(iCell)], 'Color',[0.1 0 0], 'LineWidth', 1 );
            end
            hold off;
        end

        refreshCellNumbers;


        set( mainFigure, 'Pointer', 'arrow' );
    end



%% cell selection
    function pressCell( hObj, event )
        xy = get( twoPhotonAxes, 'CurrentPoint' );
        xy=xy(end,1:2);
        if ~isempty( twoPhotonRGB ) && xy(1) >=1 && xy(2) >= 1 && xy(1) <= size(twoPhotonRGB, 2) && xy(2) <= size(twoPhotonRGB, 1)
            d = sqrt( (xy(1)-traceInfo.img_x).^2+(xy(2)-traceInfo.img_y).^2 );
            [d,j] = min(d);
            if d < traceInfo.cell_radius(j)
                newCell = traceInfo.masknum(j);
            else
                newCell = 1;   % select neuropil if clicked outside any cell
            end
            if ismember( newCell, cellList )
                cellList = setdiff( cellList, newCell );
            else
                cellList = sort([cellList, newCell]);
            end
            str = sprintf( '%d,', cellList );
            set( cellSelectionText, 'String', ['Selected cells: ' str(1:end-1)] );
        end
        if ~isempty( cellList )
            set( plotCellsButton, 'Enable', 'on' );
            set( clearCellsButton, 'Enable', 'on' );
            refreshCellNumbers;
        end
    end



    function plotSelectedCells( obj, event )
        str = sprintf( '%d,', cellList );
        plot( CircStats(key, sprintf('masknum in (%s)',str(1:end-1))) );
    end


    function clearCellSelection( obj, event )
        cellList = [];
        set( cellSelectionText, 'String', '<No cells selected>');
        set( plotCellsButton, 'Enable', 'off' );
        set( clearCellsButton, 'Enable', 'off' );
        refreshTwoPhoton;
    end


    function refreshCellNumbers
        axes( twoPhotonAxes );
        for iCell=cellList
            idx = find( iCell==traceInfo.masknum );
            ht = text( traceInfo.img_x(idx)+6, traceInfo.img_y(idx)-4, num2str( iCell), 'Color', [0 0 0], 'fontsize', 16,'fontweight','bold');
        end
    end



%% navigation
    function goToNewScan( obj, event )
        cellList = [];
        set( plotCellsButton, 'Enable', 'off' );
        set( clearCellsButton, 'Enable', 'off' );

        % lock the Scans(key) choice until the last one is fully updated
        set( scanChoice, 'Enable', 'inactive' );
        cla( twoPhotonAxes );

        % switch to the new Scans(key)
        key = scanKeys( get( scanChoice, 'value' ) );

        % list available map choices
        mapChoices = struct('code', {}, 'display', {} );
        mapChoices(end+1) = struct( 'code', 'gray'     , 'display', 'gray'                           );
        mapChoices(end+1) = struct( 'code', 'red/green', 'display', 'red/green'                      );
        set( mapChoiceMenu, 'String', {mapChoices.display} );

        % list the available parameter choices
        traceOptions = take( TraceOpts*CircStats(key),'trace_opt','brief');
        set(paramChoice,'String',{traceOptions.brief},'Value',min(max(1,get(paramChoice,'Value')),length(traceOptions)));
        if ~isempty(traceOptions)
            set(paramChoice,'Visible','on');
        else
            set(paramChoice,'Visible','off');
        end

        % refresh images
        refreshIntrinsic;
        refreshTwoPhoton;

        % clear the cell list and the comment box
        cellList = [];
        set( plotCellsButton, 'Enable', 'off' );
        set( clearCellsButton, 'Enable', 'off' );
        set( editNote, 'String', take( Scans(key), 'user_notes'));
        set( cellSelectionText, 'String', '<No cells selected>');

        % re-enable the GUI control
        set( scanChoice, 'Enable', 'on');
    end



    function printPdf( obj, event )
        filename = sprintf( './Figure_%s-%d', key.exp_date, key.scan_idx );
        print( mainFigure, '-dpdf', '-noui', filename );
        makeMovie = false;
        if makeMovie
            tpr = tpReader(getLocalPath(take(Scans(key), 'two_photon_file')));
            timeStamps = take( VisStims, 'frame_timestamps', key )/1000;
            timeStamps = timeStamps - timeStamps(1);
            makeColorMovie( tpr, sprintf( '%s_%03d', key.exp_date, key.scan_idx), find( timeStamps > 0 & timeStamps < 185 ) );
        end
    end



    function saveComment( obj, event )
        str = cellstr( get( editNote, 'String' ) );
        str = sprintf( '%s\n', str{:});
        put( Scans(key), 'user_notes', str(1:end-1) );
    end
end