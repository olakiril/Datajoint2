[cellx celly cellz] = fetchn(MaskCells(key),'img_x','img_y','img_z');
%%
f =-3;
iNeuron = 40;

cx = cellx.*exp(f);
cy = celly.*exp(f);
cz = cellz.*exp(-f);

cellDist = euDist(cx(iNeuron),cx,cy(iNeuron),cy,cz(iNeuron),cz);
[dist sortseq] = sort(cellDist);
indx = sortseq(1:key.neurons);
z = ones(size(cellx));
z(indx) = 2;
z(indx(1)) = 3;

scatter3(cellx,celly,cellz,100,z,'filled')
