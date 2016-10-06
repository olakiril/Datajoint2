%{
mov3d.RSVM (computed) # calcium trace
-> mov3d.ReduceDM
---
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
            for i = 1:7
                Mdl = fitrsvm(mapped,Par(i,:)','Standardize',true,'KFold',5,'KernelFunction','gaussian');

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
            tuple.kfoldloss = kfloss;
            tuple.rmse = rsqe;
            tuple.param_labels = labels;
            self.insert(tuple)
        end
    end
end
