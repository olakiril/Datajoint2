function makeTuples( obj, key )

[rf_snr_thr rf_p_thr bw gamma psi lambda ori_num] = ...
    fetch1(StatsOriParams(key),'rf_snr_thr', 'rf_p_thr', 'bw', 'gamma', 'psi', 'lambda', 'ori_num');

% select cells
k.mouse_id = key.mouse_id;
k.exp_date= key.exp_date;
k.scan_idx= key.scan_idx;

% get traces of all neurons of a site and fps
T = Traces(k);

if ~isempty(RFFit(k))
    T = T.*RFFit(k,['snr >' num2str(rf_snr_thr)]);
end
if ~isempty(RFStats(k))
    if rf_p_thr==0; thr = 0.05;else thr = 1;end
    T = T.*RFStats(k,['onpoff_p <' num2str(thr)]);
end
if isempty(T)
    display('Too many restrictions!, not enough RFs')
    T = Traces(k);
end

% get center of site RF
map = fetchn(RFMap*T,'onpoff_rf');
if isempty(map)
    disp('site does not have a RFMap, skipping..')
    return
end
map = map{1};
rf = fetchn(RFFit.*T,'gauss_fit');
center = mean(cell2mat(cellfun(@(x) x(1:2),rf,'UniformOutput',0)'),2);

% types of movies shown
movieTypes = unique(fetchn(StatsPresents(key),'movie_num'));

% built gabor
i = 0;
gb = repmat(gabor_fn(bw,gamma,psi,lambda,0),[1 1 ori_num]);
for iOri = 0:round(180/ori_num):180-round(180/ori_num)
    i = i+1;
    theta   = iOri/180 * pi ;
    gb(:,:,i) = gabor_fn(bw,gamma,psi,lambda,theta);
end

% read movies
if datenum(key.exp_date,'yyyy-mm-dd')>datenum('2011-12-01','yyyy-mm-dd')
    path = getLocalPath('/lab/users/philipp/stimuli/MouseMovie/MacNew/mov1l11s');
else
    path = getLocalPath('/lab/users/philipp/stimuli/MouseMovie/Mac/mov');
end
if strcmp(key.movie_type,'phase')
    type = 'phs';
else
    type = 'nat';
end
   movie = mmreader([path num2str(movieTypes(1)) '_' type '.avi']); 
 
% movie = mmreader(getLocalPath([path '/mov' num2str(movieTypes(1)) moviesEnd]));
center = round([(center(2)/size(map,1))*movie.H (center(1)/size(map,2))*movie.W]);
key.R = nan(length(movieTypes),size(gb,3),movie.NumberOfFrames);
for iMovie = 1:length(movieTypes)
    disp(['Movie: ' num2str(iMovie)])
    movie = mmreader(getLocalPath( [path num2str(movieTypes(iMovie)) '_' type '.avi'])); %#ok<TNMLP>
    for iFrame = 1:movie.NumberOfFrames
        frame = mean(read(movie,iFrame),3);
        for iFilter = 1:size(gb,3)
            [a b] = immatch(frame,gb(:,:,iFilter),center);
            key.R(iMovie,iFilter,iFrame) = dot(a(:),b(:));
        end
    end
end

% insert data
insert( obj, key );



