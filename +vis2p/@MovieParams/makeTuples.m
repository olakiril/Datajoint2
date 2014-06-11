function makeTuples( obj, key )

% get stim file
stim = fetch1(VisStims(key),'stim_file');
      
% insert
insert( obj, key );

% fill the dependant tables
makeTuples(MoviePresents,key,stim);




