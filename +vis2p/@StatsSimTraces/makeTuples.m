function makeTuples( obj, key )

% get processing parameters
options = fetch1(StatsSimTracesParams(key),'filter_type');
options = strread(options, '%s','delimiter',',');

% get movies
[mpath,m_res] = fetch1(StatsSimTracesParams(key),'movie_path','movie_resize');
cnames = {'*nat*','*phs*';'natural','phase'};
names = dir([mpath cnames{1,1}]);
movie = VideoReader([mpath names(1).name]);
frame = imresize(mean(read(movie, 1),3),m_res);
yxcenter = round(size(frame)/2);

% get simulated Rfs
params.filters = getFilters(obj,key,'movie_size',size(frame));

% get random positions
params.indx = randi(size(params.filters,1)^2 - 2*size(params.filters,1),...
    size(params.filters,3),1);

% loop through different movies shown
for iMov = 1:length(names)
    display(['Processing Movie: ' num2str(iMov) '/' num2str(length(names))])
    % loop through different classes
    for iClass = 1:size(cnames,2);names = dir([mpath cnames{1,iClass}]);
        
        movie = VideoReader([mpath names(iMov).name]); %#ok<TNMLP>
        im = mean(read(movie, 1),3);
        [~,~,params.ifra,params.ifil] = immatch(im,params.filters(:,:,1),yxcenter,'circlemask',0);
        
        % initialize
        responses = nan(size(params.filters,3),length(1:3:movie.NumberOfFrames));
        ifr = 0;
        for iFrame = 1:3:movie.NumberOfFrames; ifr = ifr+1;
            % get movie
            R = imresize(mean(read(movie, iFrame),3),m_res);
            
            % get responses
            for iopt = 1:length(options)
                R = eval([options{iopt} '(R,params)']);
            end
            % record responses
            responses(:,ifr) = R;
        end
        % fill in
        subkey = key;
        subkey.movie_type = cnames{2,iClass};
        subkey.movie_num = names(iMov).name;
        subkey.sim_traces = responses;
        
        % insert data
        insert( obj, subkey );
    end
end








