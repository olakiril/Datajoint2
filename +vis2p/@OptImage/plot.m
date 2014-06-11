function  plot(obj,varargin)

% plot(obj)
%
% Plots the intrinsic imaging data aquired with Intrinsic Imager
%
% MF 2012-09

params.sigma = 2; %sigma of the gaussian filter
params.gauss = 10; %size of gaussian window in pixels
params.vessels = 'vessels.h5'; % vessels image as backround
params.reverse = 0; % reverse the axis

params = getParams(params,varargin);

figure
set(gcf,'position',[50 200 1383 362])

keys = fetch(obj);

for ikey = 1:length(keys)
    
    set(gcf,'name',['OptMap ' keys(ikey).exp_date ' ' num2str(keys(ikey).scan_idx)])
    
    [imP vessels imA] = fetch1(OptImage(keys(ikey)),'spot_amp','vessels','spot_r2');
    
    m = cat(3,1,0,1);
    y = cat(3,1,1,0);
    b = 1-y;
    g = 1-m;
    amp = imP;
    v = imA;
    img = bsxfun(@times, amp(:,:,1), b) ...
        + bsxfun(@times, amp(:,:,2), m) ...
        + bsxfun(@times, amp(:,:,3), y) ...
        + bsxfun(@times, amp(:,:,4), g);
    
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);
    R = R - median(R(:));    R = R / quantile(abs(R(:)),0.99);
    G = G - median(G(:));    G = G / quantile(abs(G(:)),0.99);
    B = B - median(B(:));    B = B / quantile(abs(B(:)),0.99);
    
    img = cat(3,R,G,B)/2+0.5;
    img = max(0, min(1, img));
    img= rgb2hsv(img);
    
    subplot(131)
    img1 = img;
    
    % img1(:,:,2) = normalize(v);
    img1(:,:,3) = normalize(vessels);
    img1 = hsv2rgb(img1);
    imshow(img1)
    if params.reverse
        set(gca,'xdir','reverse')
    end
    subplot(132)
    img1 = img;
    
    img1(:,:,3) = normalize(v);
    img1(:,:,2) = normalize(vessels);
    
    img1 = hsv2rgb(img1);
    imshow(img1)
    if params.reverse
        set(gca,'xdir','reverse')
    end
    subplot(133)
    img1 = img;
    
    img1(:,:,2) = normalize(vessels);
    img1 = hsv2rgb(img1);
    imshow(img1)
    if params.reverse
        set(gca,'xdir','reverse')
    end
    if ikey ~= length(keys)
        pause
    end
end