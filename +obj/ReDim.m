%{
# Dimentionality reduction for Objects
-> obj.ReDimOpt
-> fuse.ScanDone
-> anatomy.Area
---
-> stimulus.Sync
mapped                      : longblob             # reduced data [time dim]
trials                      : mediumblob           # trial idx
bins                        : mediumblob           # bin idx
unit_ids                    : mediumblob           # ids of used cells
%}

classdef ReDim < dj.Computed
    %#ok<*INUSL>
    
    properties
        keySource  = (fuse.ScanDone * anatomy.Area & anatomy.AreaMembership)...
            * (obj.ReDimOpt & 'process = "yes"') ...
            & (stimulus.Sync & (stimulus.Trial &  (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"'))))
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            % setup params
            [bin,dims,gwin,reduce_func, cell_num, select_method] = fetch1(obj.ReDimOpt & key,...
                'binsize','dimensions','gauss_win','reduce_func','cell_num','select_method');
            gwin = floor(gwin/bin); % convert to time bins
            
            % get data
            [Data, Trials, Unit_ids ] = getData(self,key,bin);
            Bins = repmat([1:size(Data,2)]',1,size(Data,3));
            Trials = repmat(Trials(:)',size(Data,2),1);
            
            % select cells
            rseed = RandStream('mt19937ar','Seed',0);
            rand_idx = randperm(rseed,size(Data,1));
            Data = Data(rand_idx,:,:);
            Unit_ids = Unit_ids(rand_idx);
            switch select_method
                case 'all'
                    cell_idx = ones(size(Data,1),1);
                case 'rf'
                    % get all rfs to compute center
                    [x,y] = getDegrees(tune.DotRFMap & (fuse.ScanDone & key) & 'p_value<0.05',1);
                    mx = mean(x);
                    my = mean(y);
                    [x,y, keys] = getDegrees(tune.DotRFMap & (fuse.ScanDone & key) & 'p_value<0.05' & (anatomy.AreaMembership & key),1);
                    idx = (mx-10)<x & x<(mx+10) & (my-10)<y & y<(my+10);
                    sel_units = [keys(idx).unit_id];
                    [~,unit_idx] = intersect(Unit_ids,sel_units);
                    cell_idx = unit_idx(1:cell_num);
                case 'subsample'
                    cell_idx = 1:cell_num;
                otherwise
                    error('Cell selection method not supported')
            end
            Data = Data(cell_idx,:,:);
            
            % prefilter data
            if gwin>1
                Data = convn(Data,gausswin(gwin)','same');
            end
            
            % do it
            mappedX = feval(eval(reduce_func), Data(:,:), dims);
            
            % insert
            key.mapped = mappedX;
            key.trials = Trials(:);
            key.bins = Bins(:);
            key.unit_ids = Unit_ids(cell_idx);
            self.insert(key)
            
        end
    end
    
    methods
        
        function [Traces, Trials, Unit_ids] = getData(obj,key,bin)
            
            % get traces
            [Traces, caTimes, keys] = getAdjustedSpikes(fuse.ActivityTrace & (anatomy.AreaMembership & key),'soma'); % [time cells]
            Unit_ids = [keys.unit_id];

            % get rid of nans
            notnanidx = ~isnan(mean(Traces,2)); % faster than all
            Traces = Traces(notnanidx,:);
            caTimes = caTimes(notnanidx);
            
            % interpolate over time
            X = @(t) interp1(caTimes-caTimes(1), Traces, t, 'linear', 'extrap');  % traces indexed by time
            
            % fetch stimuli without repeats
            [flip_times, trial_idxs] = fetchn(stimulus.Trial &  ...
                (stimulus.Clip & (stimulus.Movie & 'movie_class="object3d" OR movie_class="multiobjects"')) & key,...
                'flip_times','trial_idx','ORDER BY trial_idx');
            ft_sz = cellfun(@(x) size(x,2),flip_times);
            tidx = ft_sz>=prctile(ft_sz,99);
            Trials = trial_idxs(tidx);
            flip_times = cell2mat(flip_times(tidx));
            
            % subsample traces
            Traces = permute(X(flip_times - caTimes(1)),[2 3 1]);
            if bin>0
                fps = 1/median(diff(flip_times(1,:)));
                d = max(1,round(bin/1000*fps));
                Traces = convn(Traces,ones(d,1)/d,'same');
                Traces = Traces(1:d:end,:,:);
            end
            Traces = permute(Traces,[2 1 3]); % in [cells bins trials]
        end
    end
end

