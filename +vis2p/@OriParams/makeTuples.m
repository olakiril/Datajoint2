function makeTuples( obj, key )

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% get the different stimulus conditions
uniStims = stim.params.constants.movieStat;
names = {'natural','phase'};

% loop through unique dotSizes and stimFrames 
for iMovie = 1:length(uniStims)
    key.movie_type = names{uniStims(iMovie)};
      
    % insert
    insert( obj, key );
    
    % fill the dependant tables
    makeTuples(StatsPresents,key,stim,uniStims(iMovie));
end



