function fhout = plotMask(obj,varargin)

params.bin = 1000;
params.plot = 1;
params.thr = 0.05;
params.save = false;

params = getParams(params,varargin);

global ind
global data
global pointi

import vis2p.*
if nargout>0
    params.plot = 0;
end

keys = fetch(obj);

fh = figure;
set(fh,'position',[ 416   620   824   358]);
for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
    clf
    if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        
        path = getLocalPath(fetch1(Experiments(key),'directory'));
        name = [path '/' key.exp_date ' ' num2str(key.scan_idx) ' Correlation Map'];
        file = dir([name '*']);
        
        if ~isempty(file)
            if params.plot
                imshow([path '/' file.name])
            end
        else
            % get traces
            fps = fetchn(Movies.*EphysTraces(key),'fps');
            tp = tpReader(Scans.*EphysTraces(key));
            
            traceP = fetchn(Traces(key,'masknum = 0'),'trace');
            key.masknum = fetchn(EphysTraces(key),'masknum');
            
            if key.masknum == 0
                continue
            end
            [x y r] = fetch1(MaskCells(key),'img_x','img_y','cell_radius');
            x = round(x);
            y = round(y);
            r = 2*round(r);
            if  isAligned(tp)
                traces =  tp.imChAligned{1}(y-r:y+r,x - r:x+r,:);
            else
                traces =  tp.imCh{1}(y-r:y+r,x - r:x+r,:);
            end
            
            sz = size(traces);
            
            traceC = trresize(double(permute(traces,[3 1 2])),fps,params.bin,'bin');
            traceP = trresize(traceP{1},fps,params.bin,'binsum')/params.bin*1000;
            
            [cor,p] = corr(traceP,reshape(traceC,size(traceC,1),[]));
            cor(p>params.thr) = 0;
            cor(cor<0) = 0;
            
            v = normalize(fetch1(Movies(key),'affine_green'));
            s = zeros(size(v));
            s(y-r:y+r,x - r:x+r) = normalize(reshape(cor,sz(1),sz(2)));
            h = ones(size(v));
            im = zeros(size(v,1),size(v,2),3);
            im(:,:,1) =h;
            im(:,:,2) =s;
            im(:,:,3) = v;
            subplot(2,3,[1 2 4 5])
            image(hsv2rgb(im))
            hold on
            rectangle('Position',[x-r,y-r,2*r,2*r])
            axis image
            title([key.exp_date ' scan:' num2str(key.scan_idx)])
            
            subplot(2,3,3)
            imagesc(reshape(cor,sz(1),sz(2)))
            colorbar
            hold on
            rectangle( 'Position', [[r+1,r+1]-r/2, r*[1 1]], 'Curvature', [1 1],'EdgeColor',[0 0 0]);
            
            axis image
            axis off
            title('Correlation Map');
            set(gcf,'name',[key.exp_date ' ' num2str(key.scan_idx) ' Correlation Map'])
            
            if nargout>0
                fhout = fh;
            end
            
            
            if params.save && isempty(dir([name '*']))
                print(fh,'-dpng',name)
            end
        end
        
        if length(keys) ~= 1 && ~params.save
            pause
        end
        
    elseif strcmp(fetch1(Scans(key),'scan_prog'),'AOD')
          name = fetch1( Experiments*Scans(key),'file_name' );
        dirs = dir(['M:\Mouse\' name(1:10) '*']);
        names = vertcat(dirs.name);
        [~,idir] = max(str2num(names(:,[12 13 15 16 18 19]))- str2num( name([12 13 15 16 18 19])));
        filename = ['M:\Mouse\' dirs(idir).name '\AODAcq\' name];
        fn = aodReader(filename,'Temporal');
        
        Fs = fn.Fs;
        eData = fn(:,2);
        [spikeTimes,spikeIdx] = detectSpikes(eData,Fs);
        
        %plot detected spikes
        fgh = figure(1001);
        plot(eData);
        hold on;
        plot(spikeIdx,zeros(size(spikeTimes)),'.r');
        title([key.exp_date ' ' num2str(key.scan_idx)])
        
        %%%% Identify patched Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find volumes
        u_ind = strfind(name,'_');
        
        volumename = fetchn(Scans(rmfield(key,'scan_idx'),[...
            'file_name LIKE "' name(1:u_ind(1)) ...
            '%" and scan_idx > ' num2str(key.scan_idx -2) ...
            ' and scan_idx < ' num2str(key.scan_idx + 2) ...
            ' and aim = "Stack"' ...
            ]),'file_name');
        if isempty(volumename)
            volumename = fetchn(Scans(rmfield(key,'scan_idx'),[...
                'file_name LIKE "' name(1:u_ind(1)) ...
                '%" and scan_idx > ' num2str(key.scan_idx -3) ...
                ' and scan_idx < ' num2str(key.scan_idx + 3) ...
                ' and aim = "Stack"' ...
                ]),'file_name');
        end
        
        if isempty(volumename)
            display('No volume, skipping...')
            return
        end
        volumename = ['M:\Mouse\' dirs(idir).name '\AODAcq\' volumename{1}];
        fn = aodReader(filename,'Functional');
        points = fn.coordinates;
        fn = aodReader(volumename,'Volume');
        
        % convert points to index
        pointi = nan(size(points));
        for icell = 1:size(points,1)
            pointi(icell,1) = find(roundall(points(icell,1),0.1) == roundall(fn.x,0.1));
            pointi(icell,2) = find(roundall(points(icell,2),0.1) == roundall(fn.y,0.1));
            pointi(icell,3) = find(roundall(points(icell,3),0.1) == roundall(fn.z,0.1));
        end
        
        data = fn(:,:,:,1);
        
        h = figure('NumberTitle','off','Menubar','none',...
            'Name','Find the patched cell',...
            'Position',[560 728 400 400],...
            'KeyPressFcn',@dispkeyevent);
        
        ind = 1;
        imagesc(data(:,:,ind))
        hold on
        celli = find(pointi(:,3)==ind);
        for icell = 1:length(celli)
            text(pointi(celli(icell),2),pointi(celli(icell),1),'+',...
                'HorizontalAlignment','center','VerticalAlignment','middle')
            text(pointi(celli(icell),2),pointi(celli(icell),1),num2str(celli(icell)),...
                'VerticalAlignment','top')
        end
        
        axis image
        colormap gray
        display(['cellnum: ' num2str(fetch1(EphysTraces(key),'masknum'))])
         input('Enter cell #: ','s')
        close(h)
        close(fgh)
    end
end

function dispkeyevent(~, event)
global ind
global data
global pointi

if strcmp(event.Key,'uparrow')
    if ind<size(data,3)
        ind = ind+1;
    end
elseif strcmp(event.Key,'downarrow')
    if ind>1
        ind = ind-1;
    end
end
clf
imagesc(data(:,:,ind))
axis image
hold on
celli = find(pointi(:,3)==ind);
for icell = 1:length(celli)
    text(pointi(celli(icell),2),pointi(celli(icell),1),'+',...
        'HorizontalAlignment','center','VerticalAlignment','middle')
    text(pointi(celli(icell),2),pointi(celli(icell),1),num2str(celli(icell)),...
        'VerticalAlignment','top')
end