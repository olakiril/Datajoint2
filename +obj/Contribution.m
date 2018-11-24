%{
-> obj.Dec
---
auc_total        :longblob                #area under the curve [classes reps]
auc              :longblob                #area under the curve [classes reps neurons]
%}

classdef Contribution < dj.Computed

    properties
        keySource = obj.Dec
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            steps = [0:1:100];
            Ef = [];
            [p,t,c] = fetch1(obj.Dec & key,'p','trial_info','classifier');
            
            [Data, Stims, info, Unit_ids] = getData(obj.Dec,key,[],0);
            
            dec_stims = cellfun(@unique,t.names);Stim = [];
            for istim = 1:length(dec_stims)
                stim_idx(istim) = find(strcmp(Stims,dec_stims{istim}));
                Stim = [Stim;ones(size(Data{stim_idx(istim)},2),1)*istim];
            end
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
    
     methods
       function plot(self,varargin)
            
            params.fontsize = 10;
            params.v = 1;
            params.target_cell_num = [];
            
            params = getParams(params,varargin);
             
            % get data
            [area, auc, auc_total, keys] = fetchn(self & 'brain_area <> "unknown"',...
                'brain_area','auc','auc_total');
            
            areas = unique(area);
            MI = cell(size(areas));
            colors = hsv(length(areas));
            for iarea = 1:length(areas)
                idx = (strcmp(area,areas(iarea)));
                c = cell2mat(cellfun(@(x,y) sort(squeeze(mean(mean(bsxfun(@minus,y,x),2))),'ASC'),auc(idx),auc_total(idx),'uni',0)')';
                mean(kurtosis(c,[],2))
                errorPlot(1:size(c,2),c,'errorColor',colors(iarea,:))
            end
             
            l = legend(areas);
            set(l,'box','off','location','northwest')
            ylabel('AUC gain')
            xlabel('Neurons (sorted)')
            title('Contribution')
       end
     end
    
     
end


