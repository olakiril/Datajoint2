for i = 1:size(d,3);imagesc(d(:,:,i));drawnow;pause(0.01);end

d= bsxfun(@minus,single(data),single(mean(data,3)));
d= bsxfun(@rdivide,single(d),single(mean(data,3)));
size(d)
for i = 1:size(d,3);imagesc(d(:,:,i));drawnow;end

for i = 1:size(d,3);imagesc(data(:,:,i));drawnow;pause(0.1);end

%%
data = permute(data,[3 1 2]);
im = double(data(:,:));
[u,~] = eigs(cov((im)),1);
pc = im*u*u';
im = (im - pc);
pd = reshape(im,size(data));
pc = reshape(pc,size(data));
pd = permute(pd,[2 3 1]);
pc = permute(pc,[2 3 1]);
data = permute(data,[2 3 1]);
%%
im2= bsxfun(@minus,single(data),single(mean(data,3)));
% im2= bsxfun(@rdivide,single(im2),single(mean(data,3)));

pc2= bsxfun(@minus,single(pc),single(mean(pc,3)));
% pc2= bsxfun(@rdivide,single(pc2),single(mean(pc2,3)));

pd2= bsxfun(@minus,single(pd),single(mean(pd,3)));
% pd2= bsxfun(@rdivide,single(pd2),single(mean(pd2,3)));
%%
for i = 1:size(data,3);
    subplot(3,1,1)
    imagesc(data(:,:,i));
        subplot(3,1,2)
          imagesc(pd(:,:,i),[-150 400]);
            subplot(3,1,3)
              imagesc(pc(:,:,i),[-10 400]);
    drawnow;
    pause(0.01);
end

% for i = 1:size(pc,3);imagesc(pc(:,:,i),[-10 400]);drawnow;pause(0.01);end