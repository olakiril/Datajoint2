function bw = drawCells(bw)
% lets the user outline objects in the image on the current axis

sz = axis;
sz = round([sz(4)-sz(3), sz(2)-sz(1)]);
if nargin<1 || isempty(bw)
    bw = false(sz);
end
assert(size(bw,1)==sz(1) && size(bw,2)==sz(2))

r = 3;  % dot radius
outlines = {};
undoBuffer = {bw};


printInstructions
while true
    [x,y,b] = ginput(1);
    if isempty(b)
        b = 13;  % special case for ENTER
    end
    
    switch b
        case num2cell('1':'9')
            % set brush size
            r = b-48;
        case 1
            bw(enumeratePixels) = true;
            redraw
        case 2
            bw(enumeratePixels) = false;
            redraw
        case 8
            % UNDO
            if length(undoBuffer)>=2
                bw = undoBuffer{end-1};
                undoBuffer(end-1:end) = [];
                redraw
            end
        case 32  % SPACE - toggle outlines
            if isempty(outlines)
                redraw
            else
                removeOutlines
            end
        case 27
            yes = questdlg('Do you wish to abort segmenting this image?','Finish segmentation', 'yes','no','no');
            if strcmpi('yes', yes)
                bw = [];  % destroy mask
                break
            end
        case 13
            yes = questdlg('Ready to commit?','Finish segmentation', 'yes','no','no');
            if strcmpi('yes', yes)
                break
            end
        otherwise
            disp '---'
            printInstructions
    end
end

    function printInstructions
        disp INSTRUCTIONS:
        disp 'Click to add pixels to mask'
        disp 'Shift-click to delete pixels from mask'
        disp 'Press BACKSPACE to undo'
        disp 'Press 1-9 to set brush size'
        disp 'Press SPACE to toggle outlines'
        disp 'Press ESC to discard all edits'
        disp 'Press ENTER to commit'
    end

    function idx = enumeratePixels
        idx = [];
        for xi=round(x+(-r:r))
            for yi=round(y+(-r:r))
                if (xi-x)^2 + (yi-y)^2 < r^2*0.64 && ...
                        yi>=1 && yi<=sz(1) && ...
                        xi>=1 && xi<=sz(2)
                    idx(end+1) = sub2ind(sz,yi,xi); %#ok<AGROW>
                end
            end
        end
    end


    function removeOutlines
        % delete old outlines
        for i=1:length(outlines)
            if ishandle(outlines{i})
                delete(outlines{i})
            end
        end
        outlines = {};
    end

    function redraw
        removeOutlines
        if isempty(undoBuffer)
            return
        end
        if ~all(undoBuffer{end}(:)==bw(:))
            % buffer only if modified
            undoBuffer = [undoBuffer(max(1,end-30):end) {bw}];  % append but limit buffer size
        end
        bounds = bwboundaries(bw,4);
        hold on
        for i=1:length(bounds)
            bound = bounds{i};
            outlines{i} = plot(bound(:,2), bound(:,1), 'r');
        end
        hold off
    end
end