function makeTuples( obj, key)
global ind
global data
global pointi

import vis2p.*
if  strcmp(fetch1(Scans(key),'aim'),'patching')
    if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
        %%%% Detect Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fs = fetch1(Movies(key),'ephys_fs');
        eData = fetch1(MaskTraces(key,'masknum = 0'),'ephys_trace');
        [lens,mag] = fetch1(Scans(key),'lens','mag');
        pixelPitch = 11000/512/mag/lens;
        [spikeTimes,spikeIdx] = detectSpikes(double(eData),Fs);
        
        %plot detected spikes
        fgh = figure(1001);
        plot(eData);
        hold on;
        plot(spikeIdx,zeros(size(spikeTimes)),'.r');
        title([key.exp_date ' ' num2str(key.scan_idx)])
        
        %%%% Identify patched Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = plot(MaskGroup(key));
        go = false;
        while ~go
            % Construct a questdlg with three options
            choice = questdlg('Does the patched cell have a number?', ...
                'Patch Menu', ...
                'Yes','Manualy segment','No','Yes');
            % Handle response
            switch choice
                case 'Yes'
                    prompt = {'Enter patched cell mask #:'};
                    dlg_title = 'Patched Cell #';
                    num_lines = 1;
                    def = {'0'};
                    out = inputdlg(prompt,dlg_title,num_lines,def);
                    key.masknum  = str2num(out{1});
                    close(f)
                    go = true;
                case 'Manualy segment'
                    close(f)
                    im = fetch1( Movies(key),'affine_green');
                    if isempty(im)
                        im = fetch1( Movies(key),'raw_green');
                    end
                    fh = figure;
                    imagesc(im);
                    colormap gray
                    hold on;
                    [x, y] = ginput(1);
                    radi = mean(fetchn(MaskCells(key,'masknum > 0'),'cell_radius'));
                    rectangle( 'Position', [[x,y]-radi, 2*radi*[1 1]],...
                        'Curvature', [1 1],'EdgeColor',[0 0 1]);
                    
                    choice = questdlg('Is this better?', ...
                        'Patch Menu', ...
                        'Yes','No','Yes');
                    switch choice
                        case 'Yes'
                            close(f)
                            go = true;
                            cellFinder = CellFinder( im, pixelPitch,...
                                'x', x, 'y', y, 'radii', radi,'contrast',1,...
                                'sharpness',80,'corrs',1);
                            key.masknum = max(fetchn(MaskCells(key),'masknum'))+1;
                            makeTuples(MaskTraces,key,cellFinder)
                        case 'No'
                            go = false;
                            close(fh)
                    end
                case 'No'
                    key.masknum = 0;
                    close(f)
                    go = true;
            end
        end
        
        %%%% Insert data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        key.spike_times = spikeTimes;
        insert(obj,key);
        close(fgh)
        
    elseif strcmp(fetch1(Scans(key),'scan_prog'),'AOD')
        
        %%%% Detect Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        name = fetch1( Experiments*Scans(key),'file_name' );
        dirs = dir(['M:\Mouse\' name(1:10) '*']);
        if isempty(dirs); display('no file found!');return;end
        names = vertcat(dirs.name);
        timediff = str2num(names(:,[12 13 15 16 18 19]))- str2num( name([12 13 15 16 18 19]));
        timeind = find(timediff<0);
        [~,i] = max(timediff(timeind));
        itime = timeind(i);
        filename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' name];
        
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
        volumename = ['M:\Mouse\' dirs(itime).name '\AODAcq\' volumename{1}];
        [x, y, z] = fetchn(MaskCells(key),'img_x','img_y','img_z');
        points = [x, y, z];
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
        
        key.masknum = str2num(input('Enter cell #: ','s'));
        key.spike_times = spikeTimes;
        insert(obj,key);
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


