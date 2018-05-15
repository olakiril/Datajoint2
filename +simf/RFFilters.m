%{
->simf.RFParams
filter_id                             : smallint      # 
---
filter                            : blob                          # 
%}

classdef RFFilters < dj.Manual
    methods
        function plot(self)
            [filters] = fetchn(self,'filter');
            filters = cell2mat(cellfun(@(x) x(:),filters,'uni',0)');
            figure
            imagesc(reshape(permute(reshape(filters,16,16,16,16),[1 3 2 4]),256,[]))
        end
    end
end

