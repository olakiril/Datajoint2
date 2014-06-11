function plot( obj )

keys = fetch(obj);
for ikey = 1:length(keys)
    key = keys(ikey);
    [fhMap, reversal, frev] = fetch1(StatArea(key),'fhMap','reversal','frev');
    imagesc(normalize(fhMap,2))
    hold on
    plot(reversal,1:size(fhMap,1),'.k')
    plot(frev,1:size(fhMap,1),'.y')
    title([key.exp_date ' ' num2str(key.scan_idx)])
    set(gcf,'position',[699   510   486   156])
     set(gcf,'name','Reversal')
    if ikey ~=length(keys)
        pause
        clf
    end
end