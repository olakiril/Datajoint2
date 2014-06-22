function makeTuples( obj, key )


    insert( obj, key ); 
makeTuples(vis2p.MaskTraces,key,[],[],1)

if ~strcmp(fetch1(vis2p.Scans(key),'scan_prog'),'Unirec')
    makeTuples(vis2p.MaskTracesQuality,key)
end
   
