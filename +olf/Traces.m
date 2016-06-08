%{
olf.Traces (computed) #
-> olf.Scan
-> olf.TracesOpt
masknum         : mediumint              #
---
trace                       : mediumblob                    # c) filtered traces using opt
ts=CURRENT_TIMESTAMP        : timestamp                     # c) extaction timestamp
quality=null                : float                         # m) correlation between reconstructed and raw trace
%}


classdef Traces < dj.Relvar & dj.AutoPopulate
    
    properties (Constant)
        popRel = (olf.Scan*olf.TracesOpt('process = "yes"'))
    end
    
    
    methods
        function self = Traces(varargin)
            self.restrict(varargin{:})
        end
    end
    
    methods(Access=protected)
        
        function makeTuples(obj,key)
            
            % get traces and frames
            [masknums,traces(:,:,1),traces(:,:,2)] = ...
                fetchn(olf.MaskTraces(key),'masknum','trace_ch1','trace_ch2');
            
            fps = fetch1(olf.Scan(key),'fps');
            traceOpt = fetch(olf.TracesOpt(key),'*');
            options = fetch1(olf.TracesOpt(key),'trace_computation');
            options = strread(options, '%s','delimiter',',');
            traceOpt.key = key;
            
            % filter traces
            traceTypes = size(traces,3);
            traces = double([traces{:}]);  % traces are in columns
            traces = reshape(traces,size(traces,1),[],traceTypes);
            for iopt = 1:length(options)
                [traces, qual, traceOpt]= eval([options{iopt} '(traces,fps,traceOpt)']);
            end
            
            % insert cell traces and annulus traces
            for imask=1:length(masknums)
                tuple = key;
                tuple.masknum = masknums(imask);
                tuple.trace = single(traces(:,imask));
                tuple.quality = qual(imask);
                insert( obj, tuple );
            end
        end
    end
end