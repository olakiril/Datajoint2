%{
trace_opt       : smallint unsigned      #
---
activ_func="ReSqu"                : enum('ReLu','ReSqU','biReLu','biReSqU','divNorm')       # activation function nonlinearity
%}


classdef TraceParams < dj.Lookup
    methods
        function activ_func = getActivation(self)
            switch fetch1(self,'activ_func')
                case 'ReLu'
                    activ_func = @(x) max(0,x);
                case 'ReSqU'
                    activ_func = @(x) max(0,x).^2;
                case 'biReLu'
                    activ_func = @(x) abs(x);
                case 'biReSqU'
                    activ_func = @(x) x.^2;
                case 'divNorm'
                    activ_func = @(x) bsxfun(@rdivide,max(0,x).^2,0.01+sum(max(0,x).^2));
            end
        end
    end
end

