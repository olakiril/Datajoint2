function matchCells( obj )

global markId
global stop

key = fetch(obj);

close all

mouse_strain = fetch1(Mice(key),'mouse_strain');
if strcmp(mouse_strain,'SST-Ai9'); type = 'SST';    
elseif strcmp(mouse_strain,'PV-Ai9'); type = 'PV';
elseif strcmp(mouse_strain,'VIP-Ai9'); type = 'VIP';
else disp 'Not known type, please update file!'; return;
end

plotStacks(key)

while ~stop
    pause
end

% insert data
markId = unique(markId);
for icell = 1:length(markId)
    key.masknum = markId(icell);
    update(MaskTraces(key),'mask_type',type);
end