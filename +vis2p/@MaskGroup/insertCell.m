function insertCell( obj, key,masktype)

masknum = [];
if strcmp(fetch1(Scans(key),'scan_prog'),'MPScan')
    %%%% Detect Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lens,mag] = fetch1(Scans(key),'lens','mag');
    pixelPitch = 11000/512/mag/lens;
    
    %%%% Identify patched Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = plot(MaskGroup(key));
    
    % complete key
    key = fetch(MaskGroup(key));
    go = false;
    while ~go
        % Construct a questdlg with three options
        choice = questdlg('Does the cell have a number?', ...
            'Cell Menu', ...
            'Yes','Manualy segment','Cancel','Yes');
        
        % Handle response
        switch choice
            case 'Yes'
                prompt = {'Enter patched cell mask #:'};
                dlg_title = 'Patched Cell #';
                num_lines = 1;
                def = {'0'};
                out = inputdlg(prompt,dlg_title,num_lines,def);
                masknum  = str2num(out{1});

                if ~isempty(masknum) 
                    display('Updating mask type...')
                     update(MaskTraces,fetch(MaskTraces(key,['masknum =' num2str(masknum)])),'mask_type',masktype);
                    close(f)
                    f = plot(MaskGroup(key));
                end
            case 'Manualy segment'
                if strcmp(masktype,'SST') || strcmp(masktype,'SST');  type = 'red'; 
                else type = 'green';end
                
                im = fetch1( Movies(key),['affine_' type]);
                if isempty(im)
                    im = fetch1( Movies(key),['raw_' type]);
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
                        cellFinder = CellFinder( im, pixelPitch,...
                            'x', x, 'y', y, 'radii', radi,'contrast',1,...
                            'sharpness',80,'corrs',1);
                        masknum = max(fetchn(MaskCells(key),'masknum'))+1;
                        tuple = key;
                        tuple.masknum = masknum;
                        makeTuples(MaskTraces,tuple,cellFinder,masktype)
                end
                close(fh)
                try close(f);end
                f = plot(MaskGroup(key));
            case 'Cancel'
                close(f)
                masknum = [];
                go = true;
        end
        
    end
end


