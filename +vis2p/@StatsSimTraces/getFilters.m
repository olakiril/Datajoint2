function filters = getFilters(StatsSimTraces,key,varargin) %#ok<INUSL>

% function filters = getFilters(StatsSimTraces,key,varargin)
%
% Gets the filters in [X Y filters]
%
% MF 2013-04

% set default parameters
p.movie_size = [479 719];

% update parameters if supplied
for i = 1:2:length(varargin);p.(varargin{i}) = varargin{i+1};end

% get parameters
[f_type,m_res,f_siz] = fetch1(StatsSimTracesParams(key),'filter_type','movie_resize','filter_size');

% if strfind(f_type,'sparsenet')
    Fil = load(getLocalPath('/lab/users/Manolis/Matlab/Libraries/sparsenet/A16.mat'));
    [L, M]=size(Fil.A);sz=sqrt(L);
    target_size = round(min(p.movie_size)*m_res*f_siz); % in pixels
    res = target_size/sz;
    filters = nan(size(imresize(reshape(Fil.A(:,1),sz,sz),res),1),...
        size(imresize(reshape(Fil.A(:,1),sz,sz),res),1),5);
    for i = 1:M;filters(:,:,i) = imresize(reshape(Fil.A(:,i),sz,sz),res);end    
% else
%     filters = [];
% end