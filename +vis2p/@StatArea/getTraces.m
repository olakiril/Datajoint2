function T = getTraces(obj,varargin)

% function T = getTraces(obj,varargin)
%
% gets the traces for the StatArea experiment
% [cells time trials]
%
% MF 2013-08

params.key = '';
params.area = 'v1o';
params.collapse = 0;

params = getParams(params,varargin);

skeys = fetch(StatsSites(params.key));
keys = skeys;

maskIds = fetch1(obj,params.area);
maskKey = sprintf('masknum=%d or ',maskIds);
maskKey = maskKey(1:end-3);

% collapse same movies
if length(keys)>1 && params.collapse==1
    keys = arrayfun(@(x) structfun(@(x) num2str(x),x,'uniformoutput',0),keys);
    keys = struct2cell(arrayfun(@(x) rmfield(x,'stim_idx'),keys));
    ind = [];for i = 1:size(keys,1);ind(i) = length(unique(keys(i,:)))==1;end
    if ind; keys = skeys(1);keys = rmfield(keys,'stim_idx');end
end

params = rmfield(params,'key');params = rmfield(params,'area');

T = cell(1,length(keys));
for k = 1:length(keys); key = keys(k);
    T{k} = getTraces(StatsSites(key),'key',maskKey,params);
end

% convert to array
if iscell(T) && length(T) == 1
    T = T{1};
end




