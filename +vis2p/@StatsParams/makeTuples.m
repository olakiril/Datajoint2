function makeTuples( obj, key )

% get stim file
stim = fetch1(vis2p.VisStims(key),'stim_file');

% get the different stimulus conditions
uniStims = stim.params.constants.movieStat;
names = {'natural','phase'};
if isfield(stim.params.constants,'trialTrigger')
    uniTrigger = stim.params.constants.trialTrigger;
end

% loop through unique dotSizes and stimFrames
for iMovie = 1:length(uniStims)
    key.movie_type = names{uniStims(iMovie)};
    
    if isfield(stim.params.constants,'trialTrigger')
        for iTrigger = 1:length(uniTrigger)
            
            key.trial_trigger = uniTrigger(iTrigger);
            
            % insert
            insert( obj, key );
            
            % fill the dependant tables
            makeTuples(vis2p.StatsPresents,key,stim,uniStims(iMovie),uniTrigger(iTrigger));
        end
    else
        % insert
        insert( obj, key );
        
        % fill the dependant tables
        makeTuples(vis2p.StatsPresents,key,stim,uniStims(iMovie));
    end
end



