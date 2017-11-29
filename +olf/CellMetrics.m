%{
olf.CellMetrics (computed) #
-> olf.OlfResponses
---
reliability_on                : float                    # explained variance
mean_on                       : float                    #  average cell response
responsiveness_on             : mediumblob               # odor specific responsiveness ttest
reliability_off               : float                    # explained variance
mean_off                      : float                    #  average cell response
responsiveness_off            : mediumblob               # odor specific responsiveness ttest
tuning_width                  : float                    # FWHM
peak_on                       : float                    #  peak cell on response
peak_off                      : float                    #  peak cell off response
%}


classdef CellMetrics < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel =  olf.OlfResponses
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
 
            
            [resp_on, resp_off] = fetch1(olf.OlfResponses & key,'resp_on','resp_off');
            tuple = key;
            
            % mean
            tuple.mean_on = nanmean(resp_on(:));
            tuple.mean_off = nanmean(resp_off(:));

            if all(tuple.mean_off(:)==0) || all(tuple.mean_on(:)==0);disp('No data!');return;end
            
            % reliability
            tuple.reliability_on = reliability(permute(resp_on,[3 1 2]));
            tuple.reliability_off = reliability(permute(resp_off,[3 1 2]));

            % tuning width
            resp = sort(nanmean(resp_on,2),'descend');
            [xData, yData] = prepareCurveData( [],  double(resp-min([0 min(resp(:))])));
            ft = fittype( 'exp1' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [0 0];
            [fitresult, ~] = fit( xData, yData, ft, opts );
            tuple.tuning_width = -log(fitresult.a/(yData(1)/2))/fitresult.b;

            if isinf(tuple.tuning_width) || isnan(tuple.tuning_width);
                print('Data problem!');return;end
            % responsivess
            tuple.responsiveness_on = ttest(resp_on',0);
            tuple.responsiveness_off = ttest(resp_off',0);
            
            % peak response
            tuple.peak_on = max(nanmean(resp_on,2));
            tuple.peak_off = max(nanmean(resp_off,2));

            % insert
            insert( obj, tuple );
            
        end
    end
end