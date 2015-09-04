function [traces,out,traceOpt] = pcMinus2(traces,fps,traceOpt)

for i = 1:size(traces,3);
    traces = double(traces);

    % if isfield(traceOpt,'highPass')
    %     if traceOpt.highPass>0
    %         k = hamming(round(fps/traceOpt.highPass)*2+1);
    %         k = k/sum(k);
    %         % make sure everything is positive
    %         traces = traces + abs(min(traces(:)))+eps;
    %         traces = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass
    %         
    %         % whatever is nan get rid of it by replacing it with 0
    %         traces(isnan(traces)) = 0;
    %         
    %         % disable further filtering
    %         traceOpt = rmfield(traceOpt,'highPass');
    %     else
    %         traces = bsxfun(@rdivide,traces,mean(traces))-1;  % dF/F where F is mean
    %     end
    % end

    [c, p] = princomp(traces(:,:,i));
    traces(:,:,i) = p(:,2:end)*c(:,2:end)';
end

if nargout>1
    out = [];
end


