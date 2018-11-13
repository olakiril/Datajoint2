%{
-> obj.Dec
---
auc_total        :double                #area under the curve [classes reps]
auc              :double                #area under the curve [classes reps neurons]
%}

classdef Contribution < dj.Computed

    properties
        keySource = obj.Dec
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            steps = [0:1:100];
            Ef = [];
            [p,t,c] = fetch1(obj.Dec & k,'p','trial_info','classifier');
            
            [Data, Stims, info, Unit_ids] = getData(obj.Dec,k,[],0);
            Stim = reshape(repmat([1 2 3 4],size(Data{stim_idx(1)},2),1),[],1);
            Data = cell2mat(Data(stim_idx));
            
            for iclass = 1:length(c{1})
                for irep = 1:size(p{1},1)
                    [~,~,cell_idxB] = intersect(t.units{irep}{1},Unit_ids,'stable');
                    
                    SVM = c{irep}{iclass}{1};
                    dat = Data(cell_idxB,:);
                    AUC = [];
                    for ineuron = 1:size(SVM,2)
                        idx = true(size(SVM,2),1);
                        idx(ineuron) = false;
                        svm = SVM(:,idx);
                        Dist = @(X) dot((X - svm(3,:)')./svm(4,:)' ,repmat(svm(1,:)',1,size(X,2))) + svm(2,1);
                        d  = Dist(dat(idx,:));
                        
                        X = prctile(d,steps);Y = [];
                        for i = 1:length(X)
                            Y(i) = mean(Stim'==iclass & d>X(i));
                        end
                        AUC(ineuron) = sum(Y/mean(Stim==iclass))/length(Y);
                    end
                    
                    svm = SVM;
                    Dist = @(X) dot((X - svm(3,:)')./svm(4,:)' ,repmat(svm(1,:)',1,size(X,2))) + svm(2,1);
                    d  = Dist(dat);
                    
                    X = prctile(d,steps);Y = [];
                    for i = 1:length(X)
                        Y(i) = mean(Stim'==iclass & d>X(i));
                    end
                    AUC_all(iclass,irep) = sum(Y/mean(Stim==iclass))/length(Y);
                    Ef(iclass,irep,:) = AUC;
                end
            end
            
            % insert
            tuple = key;
            tuple.auc_total = AUC_all;
            tuple.auc = Ef;
            disp 'Inserting ...'
            insert(self,tuple)

        end
    end
end


