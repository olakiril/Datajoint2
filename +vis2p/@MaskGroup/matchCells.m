function matchCells(obj, AODKey1,AODKey2,NewMasknum)

% matches the cellids between different AOD volumes

global markId1
global markId2
global stop

close all

compareAODVolumes(vis2p.MaskGroup,AODKey1,AODKey2)

while ~stop
    pause(1)
end

% insert data
assert( ~isempty(markId1) && ~isempty(markId2),'not enough cells selected!')    

% get unique masknumbers
[~,idx1] = unique(markId1);
[~,idx2] = unique(markId2);
markId1 = markId1(intersect(idx1,idx2));
markId2 = markId2(intersect(idx1,idx2));

% get all masknumbers from volume to be updated
disp 'Fetching tables to be updated...'
MaskTracesData = fetch(vis2p.MaskTraces(AODKey2),'*');
MaskCellsData = fetch(vis2p.MaskCells(AODKey2),'*');
MaskTracesQData = fetch(vis2p.MaskTracesQuality(AODKey2),'*');
masks = [MaskTracesData.masknum];
delQuick(vis2p.MaskTraces(AODKey2));
delQuick(vis2p.MaskCells(AODKey2));
delQuick(vis2p.MaskTracesQuality(AODKey2));
disp 'done!'

% disable key constrains
disp 'disabling contrains'
query(dj.conn, 'set foreign_key_checks=0')
disp done!

% loop through all the masks and update ids
disp 'updating masknumbers'
for icell = 1:length(masks)

    mask = masks(icell);
    if any(mask==markId2) % set matching mask number
        masknum = markId1(mask==markId2);
    else % set unique number
        NewMasknum = NewMasknum+1;
        masknum = NewMasknum;
    end
    key = MaskTracesData(icell);
    key.masknum = masknum;
    insert(vis2p.MaskTraces,key);
    key = MaskCellsData(icell);
    key.masknum = masknum;
    insert(vis2p.MaskCells,key);
    key = MaskTracesQData(icell);
    key.masknum = masknum;
    insert(vis2p.MaskTracesQuality,key);
end
disp done!

% enable all key constrains
disp 'enabling contrains'
query(dj.conn, 'set foreign_key_checks=1')
disp done!
