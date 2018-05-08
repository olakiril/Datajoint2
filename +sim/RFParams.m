%{
rf_opt       : smallint unsigned      #
---
x_size                            : smallint                            # x filter size
y_size                            : smallint                            # y filter size
type                              : enum('pixels','ICA','multiICA')        # filter type
brief="fill out"                  : varchar(127)                        # short description, to be displayed in menus
process="yes"                     : enum('no','yes')                    # do it or not
%}


classdef RFParams < dj.Lookup
    methods
        function createFilters(self)
            [type,x_sz,y_sz] = fetch1(self,'type','x_size','y_size');
            switch type
                case 'pixels'
                    filters = reshape(diag(ones(x_sz*y_sz,1)),y_sz,x_sz,x_sz*y_sz);
                case 'ICA'
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                case 'multiICA'
                    filt = load(getLocalPath('/stor02/Stimuli/ICAFilters/A16.mat'));
                    filters = reshape(filt.A,16,16,256);
                    c = corr(abs(filt.A));
                    c(logical(diag(ones(size(c,1),1))))=nan;
                    for ifilter = 1:size(filters,3)
                        [~,idx] = nanmax(c(ifilter,:));
                        filters(:,:,ifilter) = filters(:,:,ifilter) + filters(:,:,idx);
                    end
            end
            
            % insert filters
            key = fetch(self);
            for ifilter = 1:size(filters,3)
                key.rf_id = ifilter;
                key.filter = filters(:,:,ifilter);
                insert(sim.RFFilters, key)
            end
        end
    end
end

