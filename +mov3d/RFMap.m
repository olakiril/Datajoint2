%{
mov3d.RFMap (computed) # population RF
-> rf.Sync
---
map                    : longblob        # standard deviation from center of population RF
rf_idx                 : longblob        # index of subbins with stimulus in pop RF {trials}(subbins)               
rf_trials              : longblob        # trial index (trial indexes)
%}

classdef RFMap < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = (rf.Scan & 'cortical_area = "V1"') * (pre.SpikeInference & 'spike_inference = 2')...
            * (pre.SegmentMethod & pre.Trace) * ...
            (rf.Sync & (psy.MovieInfo & 'movie_class="object3d"') & psy.MovingNoise) & ...
            pre.Spikes & monet.Fit
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            % binsize = 
            
            [xloc, yloc, r, c] = fetchn(monet.Fit & key,'x','y','radius','contrast');
            [width, height] = fetchn(monet.RF & key, 'degrees_x', 'degrees_y');
            m = fitgmdist([xloc,yloc],1);
            mu = m.mu;
            C = m.Sigma;
            [x,y] = meshgrid(-125:124,-125:124);

            X=[x(:) y(:)];

            X = bsxfun(@minus, X, mu);
            d = sum((X /C) .* X, 2);

           % stimulus_trial_xy_position
           % stimulus_bin_xy_position inside circle
 
 

            self.insert(tuple)
            
        end
    end
    
end
