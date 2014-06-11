function preview(key)

k = fetch(key);

for ikey = 1:length(k)
    if ikey ~= 1
            pause
    end
    clf
    kk = k(ikey)%#ok<NOPRT>

    [x y z id] = fetchn(MaskCells(kk),'img_x','img_y','img_z','masknum');


    scatter3(x,y,z,1,'filled')
    hold on
    for iid = 1:length(id)
         text(x,y,z,num2str(id))
    end
    set(gcf,'Color',[1 1 1])
        drawnow

        

end