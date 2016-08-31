function [J] = nnclassSVfast(traces,varargin)

% function [J rJ] = nnclass(traces)
%
% performs a nearest neighbor classification
% traces: [cells classes trials]
%
% MF 2011-08-25

params.repetitions = 1;
params.cells = size(traces,1);
params.frames = 0;
params.trials = 0;

params = getParams(params,varargin);

if params.frames; y=size(traces,2);else y=1;end

if params.cells == size(traces,1)
    params.repetitions = 1;
end

nclasses = 2;

J = cell(size(traces,3),1);
parfor iTrial = 1:size(traces,3)
    ind = true(size(traces,3),1);
    ind(iTrial) = false;
    r = traces(:,:,ind);
    group = repmat((1:nclasses)',1,size(r,3));
    SVMStruct = fitcsvm(r(:,:)',group(:));
    p =[];
    for iClass = 1:size(traces,2)
        indx = predict(SVMStruct,traces(:,iClass,iTrial)');
        p(iClass) = indx == iClass;
    end
    J{iTrial} = p;
end


