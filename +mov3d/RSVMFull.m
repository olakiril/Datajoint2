%{
mov3d.RSVMFull (computed) # calcium trace
-> mov3d.Decode
---
params                      : longblob
predict                     : longblob
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
            predict = nan(size(Par));
            for i = 1:7
                Mdl = fitrsvm(Data,Par(i,:)','Standardize',true,'KFold',5,'KernelFunction','linear');
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
            
            sessions = fetch(experiment.Session & (obj & 'dec_opt = 32'));
            areas =  fetchn(map.Area,'area');
            
            RMSE = nan(length(areas),length(sessions));
            
            for isession = 1:length(sessions)
                for iarea = 1:length(areas)
                    keys = fetch(obj & sessions(isession) & (experiment.Scan & ['brain_area="' areas{iarea} '"']));
                    if isempty(keys);continue;end
                    rmse = [];
                    for ikey = 1:length(keys)
                        tuple = keys(ikey);
                        [RM, labels] = (fetch1(obj & tuple,'rmse','param_labels'));
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
            set(gcf,'name',['Normalized RMSE ' labels{iparam}])
            c.Ticks = 0:0.1:1;
            c.TickLabels = roundall(mn_nrmse: (mx_nrmse - mn_nrmse)/5 :mx_nrmse,0.01);
            ylabel(c, 'Normalized RMSE to V1')
            
            
        end
        
        function plotBars(obj,iparam)
            [rmse,labels,area]= fetchn(obj*experiment.Scan,'rmse','param_labels','brain_area');
            
            data = cellfun(@(x) x(:,iparam),rmse);
            
            idx = ~isnan(data);
            
            data = data(idx);
            area = area(idx);
            
            uniarea = unique(area);
            
            Data = cell(length(uniarea),1);
            Area = cell(length(uniarea),1);
            for i = 1:length(uniarea)-1
                
                Area{i} = uniarea{i};
                Data{i} = data(strcmp(area,Area{i}));
                
                
            end
            figure
            barfun(Data(1:end-1)')
            
            legend(Area(1:end-1)')
            title(labels{1}{iparam})
            
            
            
        end
        
        function plotRaw(obj,varargin)

            params.
            
            params = getParams(params,varargin);

            
            [params_in,predict,labels,areas]= fetchn(obj*experiment.Scan,'params','predict','param_labels','brain_area');
            
            if nargin<2
                params = param;
            else
                params = 1:size(predict{1},2);
            end
            
            uniarea = unique(areas);
            
            for iparam = params
                for iarea = 1:length(uniarea)-1
                    fig = figure;
                    hold on
                    idx = strcmp(areas,uniarea{iarea});
                    
                    
                    dataPr = cellfun(@(x) x(iparam,:),predict(idx),'uni',0);
                    dataPa = cellfun(@(x) x(iparam,:),params_in(idx),'uni',0);
                    
                    colors = hsv(length(dataPr));
                    for isite = 1:length(dataPr)
%                        plot(dataPa{isite},dataPr{isite},'.','color',colors(isite,:)); 
%                         
                         regressPlot(dataPa{isite},dataPr{isite},'globalcolor',...
                             colors(isite,:),'figure',fig,...
                             'MarkerSize',1,'MarkerType','.'); 

                    end
                    
                    name = sprintf('%s, %s',uniarea{iarea},labels{1}{iparam});
                    title(name)
                    set(gcf,'name',name)

                end
            end
            
        end
    end
end