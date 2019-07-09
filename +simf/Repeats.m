%{
# Movie repeats evaluation
-> simf.RepeatsOpt
-> simf.Activity
%}

classdef Repeats < dj.Computed
    
    
    properties
        keySource = simf.Activity * (simf.RepeatsOpt & 'process = "yes"')
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            % get parameters & data
            disp 'fetching data ...'
            [traces,keys] = fetchn(simf.ActivityTrace & key,'trace');
            Data = cell2mat(traces);
            [method, noise, sample_sz] = fetch1(simf.RepeatsOpt & key,'method','noise','sample_sz');
            
            noise_func = eval(noise);
            
            nData = nan(size(Data,1),size(Data,2),sample_sz);
            for i = 1:sample_sz
                nData(:,:,i) = noise_func(Data);
            end
            Data = permute(nData,[2 3 1]);

            % insert
            tuple = key;
            self.insert(tuple)
            
            % initialize
            disp 'computing reliabity ...'
            import_keys = cell(size(Data,3),1);
            rl = cell(size(Data,3),1);
            for icell = 1:size(Data,3)
                
                trace = Data(:,:,icell); %[trials time]
                switch method
                    case 'explainedVar'
                        rl{icell} = var(nanmean(trace,2))/nanvar(trace(:));
                    case 'corr'
                        traceZ = zscore(trace);
                        c = corr(traceZ);
                        sel_idx = logical(tril(ones(size(traceZ,2)),-1));
                        rl{icell} = nanmean(c(sel_idx));
                    case 'oracle'
                        traceZ = zscore(trace);
                        
                        c = nan(size(traceZ,2),1);
                        for irep = 1:size(traceZ,2)
                            idx = true(size(traceZ,2),1);
                            idx(irep) = false;
                            c(irep) = corr(traceZ(:,irep),mean(traceZ(:,idx),2));
                            
                        end
                        rl{icell} = nanmean(c);
                end
                
                % clear tuple
                tuple = keys(icell);
                tuple.rep_opt = key.rep_opt;
                tuple.r = rl{icell};                
                import_keys{icell} = tuple;
            end
            
            % insert
            disp 'Inserting ...'
            for tuple = cell2mat(import_keys')
                insert(simf.RepeatsUnit,tuple)
            end
        end
    end
    
end
