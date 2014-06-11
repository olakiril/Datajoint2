function makeTuples( obj, key )


    insert( obj, key ); 
makeTuples(MaskTraces,key,[],[],1)

if ~strcmp(fetch1(Scans(key),'scan_prog'),'Unirec')
    makeTuples(MaskTracesQuality,key)
end
   
