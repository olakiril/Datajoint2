function [traces,out,traceOpt] = globalPCminus(traces2,fps,traceOpt)

[lens,mag] = fetch1(Scans(traceOpt.key),'lens','mag');
pixelPitch = 11000/512/mag/lens;
tp = tpReader(Scans(traceOpt.key));
imChan = getAlignedImagingChannel(tp,1);

imChan = imChan(:,:,:);
imChan = imresize(imChan,0.5);
imChan = permute(imChan,[3 1 2]);
dsz = size(imChan);
imChan = double(imChan(:,:));

if isfield(traceOpt,'highPass')
    if traceOpt.highPass>0
        k = hamming(round(fps/traceOpt.highPass)*2+1);
        k = k/sum(k);
        % make sure everything is positive
        imChan = imChan + abs(min(imChan(:)))+eps;
        imChan = imChan./convmirr(imChan,k)-1;  %  dF/F where F is low pass
        
        % whatever is nan get rid of it by replacing it with 0
        imChan(isnan(imChan)) = 0;
        
        % disable further filtering
        traceOpt = rmfield(traceOpt,'highPass');
    else
        imChan = bsxfun(@rdivide,imChan,mean(imChan))-1;  % dF/F where F is mean
    end
end

[u,~] = eigs(cov((imChan)),1);
imChan = (imChan - imChan*u*u');
imChan = reshape(imChan,dsz);
imChan = permute(imChan,[2 3 1]);
[x, y, radi, cntr, sharp] = fetchn(MaskCells(traceOpt.key,...
    'masknum <> -1 and masknum <> -3 and masknum <> 0'),...
    'img_x','img_y','cell_radius','green_contrast','sharpness');
cellFinder = CellFinder( mean(imChan,3), pixelPitch,...
    'x', x/2, 'y', y/2, 'radii', radi/2,'contrast',cntr,...
    'sharpness',sharp);

[traces, ntrace] = getTraces( cellFinder, imChan );
traces = [ntrace traces];


if nargout>1
    out = [];
end

%%
%%
% traces2 = traces2(:,:,1);
% traces2 = traces2 + abs(min(traces2(:)))+eps;
% traces2 = traces2./convmirr(traces2,k)-1;  %  dF/F where F is low pass
% 
% %%
% figure
% hold on
% plot(bsxfun(@plus,traces,1:size(traces,2)),'k')
% plot(bsxfun(@plus,traces2,(1:size(traces,2)))-0.4,'r')
% imChan3 = reshape(imChan3,dsz);
% imChan3 = permute(imChan3,[2 3 1]);
% %%
% 
%  for i = 1:size(imChan,3);
%       clf;
%      
%      subplot(3,1,1)
% 
%      imagesc(imChan(:,:,i),[-0.75 1]);
%       axis off
%      title(num2str(mod(i/fps,4)))
%      subplot(3,1,2)
%      imagesc(imChan2(:,:,i),[-0.75 1]);
%      axis off
%       subplot(3,1,3)
%      imagesc((imChan3(:,:,i)),[-0.75 1]);
% %      colorbar
% axis off
%      drawnow;
%      pause(0.01);
%  end