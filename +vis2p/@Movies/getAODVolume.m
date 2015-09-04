function volume = getAODVolume(self,key1)

import vis2p.*

name = fetch1( Experiments*Scans(key1),'file_name' );
if ~strcmp(fetch1(Scans(key1),'aim'),'stack')
    display('Detecting volume..')
    u_ind = strfind(name,'_');
    
    volumename = fetchn(Scans(rmfield(key1,'scan_idx'),[...
        'file_name LIKE "' name(1:u_ind(1)) ...
        '%" and scan_idx > ' num2str(key1.scan_idx -2) ...
        ' and scan_idx < ' num2str(key1.scan_idx + 2) ...
        ' and aim = "Stack"' ...
        ]),'file_name');
    if isempty(volumename)
        volumename = fetchn(Scans(rmfield(key1,'scan_idx'),[...
            'file_name LIKE "' name(1:u_ind(1)) ...
            '%" and scan_idx > ' num2str(key1.scan_idx -3) ...
            ' and scan_idx < ' num2str(key1.scan_idx + 3) ...
            ' and aim = "Stack"' ...
            ]),'file_name');
    end
    
    if isempty(volumename)
        display('No volume, skipping...')
        return
    else
        display(['Using: ' volumename{1}])
    end
    volumename = volumename{1};
else
    volumename = name;
end

% select an earlier nearby directory and load the volume
dirs = dir(['M:\Mouse\'  volumename(1:10) '*']);
dif = bsxfun(@minus,datenum(vertcat(dirs.name),'yyyy-mm-dd_HH-MM-SS'), datenum(volumename(1:19),'yyyy-mm-dd_HH-MM-SS'));
dirs = dirs(dif<0);
[~,idx] = max(dif(dif<0));
volumename = ['M:\Mouse\' dirs(idx).name '\AODAcq\' volumename];
[x, y, z, Ids] = fetchn(MaskCells(key1),'img_x','img_y','img_z','masknum');
points = [x, y, z];
fn = aodReader(volumename,'Volume','ch1');
aod_step_size = (fn.z(end) - fn.z(1))/length(fn.z);

% convert points to index
pointi = nan(size(points));
for icell = 1:size(points,1)
    [~,pointi(icell,1)] = min(abs(points(icell,1)-fn.x));
    [~,pointi(icell,2)] = min(abs(points(icell,2)- fn.y));
    [~,pointi(icell,3)] = min(abs(points(icell,3)-fn.z));
end

volume = fn(:,:,:,:);