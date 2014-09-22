function plot(obj,key)

[m,rz] = fetch1(RFMapRaw(key),'onpoff_rf','rz');
m = double(imresize(m,1/rz));

for i = 1:length(m);
    [~,X] = max( squeeze(mean(m{i},3)),[],3);
    [~,Y] = max( squeeze(mean(m{i},4)),[],3);
    clf;
    set(gcf,'position',[680   756   660   342])
    subplot(2,4,1:3)
    imagesc(X);
    title('Azimuth')
%     axis off
    c1 = colorbar;
    set(c1,'location','south','position',[0.7621    0.6374    0.1904    0.2219])
    set(c1,'xtick',[])
    subplot(2,4,5:7)
    imagesc(Y);
    title('Elevation')
%     axis off
    c2 = colorbar;
    set(c2,'location','east','position',[ 0.7621    0.1637    0.1904    0.2219])
    set(c2,'ytick',[],'ydir','reverse')
     set(gcf,'name','RFMap')
    if length(m)~=i
        pause;
    end
end