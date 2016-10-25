%{
mov3d.RSVM (computed) # calcium trace
-> mov3d.ReduceDM
---
params                      : longblob
predict                     : longblob
kfoldloss                   : longblob             # generalization error 
rmse                        : longblob             # root mean square error
param_labels                : longblob             # param label
%}

classdef RSVM < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = mov3d.ReduceDM
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            labels = {'Xpos','Ypos','Scale','Xrot','Yrot','Speed','Object'};

            [~,~, Params] = getData(mov3d.ReduceDM,key);
            mapped = fetch1(mov3d.ReduceDM & key,'mapped');
            Params(7,1,:) = ones(size(Params,3),1);
            Params(7,2,:) = 2*ones(size(Params,3),1);
            Par = Params(:,:);

            kfloss = [];
            rsqe = [];
            predict = nan(size(Par));
            for i = 1:7
                Mdl = fitrsvm(mapped,Par(i,:)','Standardize',true,'KFold',5,'KernelFunction','gaussian');
                predict(i,:) = Mdl.kfoldPredict;
                kfloss(i) = kfoldLoss(Mdl);

                vp = Mdl.kfoldPredict;
                vp = vp- mean(vp);
                vp = vp/std(vp);

                v = Par(i,:)';
                v = v- mean(v);
                v = v/std(v);

                rsqe(i) = sqrt(mean((vp-v).^2));
            end
            
            % insert
            tuple.params = Par;
            tuple.predict = predict;
            tuple.kfoldloss = kfloss;
            tuple.rmse = rsqe;
            tuple.param_labels = labels;
            self.insert(tuple)
        end
    end
    methods
        function plotMasks(obj,iparam)

            sessions = fetch(experiment.Session & (obj & 'red_opt = 1'));
            areas =  fetchn(map.Area,'area');
            
            RMSE = nan(length(areas),length(sessions));
            
            for isession = 1:length(sessions)
                for iarea = 1:length(areas)
                    keys = fetch(obj & sessions(isession) & (experiment.Scan & ['brain_area="' areas{iarea} '"']));
                    if isempty(keys);continue;end
                    rmse = [];
                    for ikey = 1:length(keys)
                        tuple = keys(ikey);
                        RM = (fetch1(obj & tuple,'rmse'));
                        rmse(ikey) = RM(iparam);
                    end
                    RMSE(iarea,isession) = mean(rmse);
                end
            end
            nrmse = nanmean(bsxfun(@rdivide,RMSE,RMSE(1,:)),2);
            mx_nrmse = max(nrmse);
            mn_nrmse = min(nrmse);
            
            figure
            colors = parula(100);
            plotMask(map.Masks)
            colormap parula
            c = colorbar;
            for iarea = 1:length(areas)
                if isnan(nrmse(iarea));continue;end
                idx = ceil(((nrmse(iarea) - mn_nrmse) /(mx_nrmse  - mn_nrmse))*100);
                if idx ==0;idx =1;end
                plotMask(map.Masks & ['area="' areas{iarea} '"'],colors(idx,:))
            end
            set(gcf,'name','Normalized RMSE')
            c.Ticks = 0:0.2:1;
            c.TickLabels = roundall(mn_nrmse: (mx_nrmse - mn_nrmse)/5 :mx_nrmse,0.1);
            ylabel(c, 'Normalized RMSE to V1')
            
            
        end
    end
end
