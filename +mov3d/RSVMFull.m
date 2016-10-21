%{
mov3d.RSVMFull (computed) # calcium trace
-> mov3d.Decode
---
kfoldloss                   : longblob             # generalization error 
rmse                        : longblob             # root mean square error
param_labels                : longblob             # param label
%}

classdef RSVMFull < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel  = mov3d.Decode & 'dec_opt = 32'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            tuple = key;
            
            labels = {'Xpos','Ypos','Scale','Xrot','Yrot','Speed','Object'};

            [Data,~,~,~,~,Params] = getData(mov3d.Decode,key);
            Params(7,1,:) = ones(size(Params,3),1);
            Params(7,2,:) = 2*ones(size(Params,3),1);
            Par = Params(:,:);

            Data = Data(:,:)';
            
            kfloss = [];
            rsqe = [];
            for i = 1:7
                Mdl = fitrsvm(Data,Par(i,:)','Standardize',true,'KFold',5,'KernelFunction','linear');

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
