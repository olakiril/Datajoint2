function plotGlobalMap(obj,key)

[m,rz] = fetchn(vis2p.RFMapRaw(key),'onpoff_rf','rz');

for i = 1:length(m);
    im = double(imresize(m{i},1/rz));

   map = squeeze(mean(mean(im,1),2));
   imagesc(map)
   axis image
   axis off
    if length(m)~=i
        pause;
    end
end