function makeTuples( obj, key )

[gwin, glim ] = fetch1(StatAreaParams(key),'gwin','glim');

rfMap = fetch1(RFMapRaw(key),'onpoff_rf');
[~,hMap] = max( squeeze(mean(rfMap,3)),[],3);
fhMap = convn(hMap,gausswin(gwin)*gausswin(gwin)','same');
[~,reversal] = max(fhMap,[],2);
fitresult = createFit(reversal);
frev = mean(predint(fitresult,1:size(hMap,1)),2);
[xCell, yCell, mnums] = fetchn(MaskCells(key,'masknum>0'),'img_x','img_y','masknum');
xlim = frev(round(yCell));

key.v1o = mnums(xCell <= xlim - glim                       );
key.v1i = mnums(xCell <  xlim        & xCell > xlim - glim );
key.v2o = mnums(xCell >= xlim + glim                       );
key.v2i = mnums(xCell >  xlim        & xCell < xlim + glim );
key.fhMap = fhMap;
key.reversal = reversal;
key.frev = frev;

% insert data
insert( obj, key );