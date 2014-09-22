function matchCells( obj,MpScanKey )

global markId
global stop

key = fetch(obj);

close all

if nargin<2
    MpScanKey = [];
end

mouse_strain = fetch1(vis2p.Mice(key),'mouse_strain');
if strcmp(mouse_strain,'SST-Ai9'); type = 'SST';    
elseif strcmp(mouse_strain,'PV-Ai9'); type = 'PV';
elseif strcmp(mouse_strain,'VIP-Ai9'); type = 'VIP';
else type = 'red';
end

if strcmp(fetch1(vis2p.Scans.*obj,'scan_prog'),'AOD')
    compareVolumes(vis2p.MaskGroup,key,MpScanKey)
else
    compareChannels(vis2p.MaskGroup,key)
end

while ~stop
    pause(1)
end

% insert data
if ~isempty(markId)
    markId = unique(markId);
    for icell = 1:length(markId)
        key.masknum = markId(icell);
        update(vis2p.MaskTraces(key),'mask_type',type);
    end
end