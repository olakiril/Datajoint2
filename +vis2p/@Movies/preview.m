function preview(key)

k = fetch(key);

for ikey = 1:length(k)
    
    kk = k(ikey)%#ok<NOPRT>
    s = fetch1(Movies(kk),'raw_green');
    figure
    for iZ = 1:size(s,3)
        imagesc(s(:,:,iZ))
        colormap gray
        drawnow
        pause
    end
end