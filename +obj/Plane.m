%{
-> obj.Dec
---
weights               :longblob                # regression weights
angle_in              :longblob                # angle between class pair [classepair reps]
angle_diag            :longblob                # angle from diagonal [classes reps]
%}

classdef Plane < dj.Computed

    properties
        keySource = obj.Dec
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            [p,t,c] = fetch1(obj.Dec & key,'p','trial_info','classifier');

            CC = cell2mat(cellfun(@(x) cell2mat(cellfun(@(xx) permute(xx{1}(1,:),[3 1 2]),x,'uni',0)'),c,'uni',0));

            for irep = 1:size(p{1},1)

                for iclass = 1:length(c{1})
                    ANGD(iclass,irep) = min(abs(acosd(dot(c{irep}{iclass}{1}(1,:),ones(size(c{irep}{iclass}{1}(1,:))))...
                        /(norm(c{irep}{iclass}{1}(1,:)))) - [0 180]));
                end
                cmb = nchoosek([1:iclass],2);
                for iclass = 1:size(cmb,1)
                    ANG(iclass,irep) = min(abs(acosd(dot(c{irep}{cmb(iclass,1)}{1}(1,:),c{irep}{cmb(iclass,2)}{1}(1,:))...
                        /(norm(c{irep}{cmb(iclass,1)}{1}(1,:))*norm(c{irep}{cmb(iclass,2)}{1}(1,:)))) - [0 180]));
                end
            end
    
            %%
            % insert
            tuple = key;
            tuple.weights = CC;
            tuple.angle_in = ANG;
            tuple.angle_diag = ANGD;
            disp 'Inserting ...'
            insert(self,tuple)

        end
        
        
    end
    
    methods
       function plot(self,type,varargin)
            
            params.fontsize = 10;
            params.v = 1;
            params.target_cell_num = [];
            
            params = getParams(params,varargin);
             
            if nargin<2 || isempty(type)
                type = 'angle_in';
            end
            
            % get data
            [area, values, keys] = fetchn(self & 'brain_area <> "unknown"',...
                'brain_area',type);
            
            areas = unique(area);
            MI = cell(size(areas));
            for iarea = 1:length(areas)
                idx = (strcmp(area,areas(iarea)));
%                   MI{iarea} = cell2mat(cellfun(@(x) x(:),values(idx),'uni',0));
               MI{iarea} =cellfun(@(x) nanmean(x(:)),values(idx));
            end
             
            boxfun(MI,'sig',1)
            set(gca,'xticklabel',areas)
       end
       
        function plotWeights(self,varargin)
            
            params.fontsize = 10;
            params.v = 1;
            params.target_cell_num = [];
            
            params = getParams(params,varargin);
             
            % get data
            [area, weights, keys] = fetchn(self & 'brain_area <> "unknown"',...
                'brain_area','weights');
            
            areas = unique(area);
            MI = cell(size(areas));
            colors = hsv(length(areas));
            for iarea = 1:length(areas)
                idx = (strcmp(area,areas(iarea)));
                c = cell2mat(cellfun(@(x) squeeze(mean(sort((mean(x,2)),3))),weights(idx),'uni',0)')';
                mean(kurtosis(c,[],2))
                errorPlot(1:size(c,2),c,'errorColor',colors(iarea,:))
            end
             
            l = legend(areas);
            set(l,'box','off','location','northwest')
            ylabel('Coefficient')
            xlabel('Neurons (sorted)')
            title('Coefficient')
        end
    end
end


