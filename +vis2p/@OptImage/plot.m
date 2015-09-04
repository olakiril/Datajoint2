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

import vis2p.*

keys = fetch(obj);

for ikey = 1:length(keys)
    
    set(gcf,'name',['OptMap ' keys(ikey).exp_date ' ' num2str(keys(ikey).scan_idx)])
    
    [imP,vessels,imA,stim] = fetch1(OptImage(keys(ikey)),'spot_amp','vessels','spot_r2','stim_loc');
    map = [];
    for istim = 1:2
        un_imP = [];
        uniq = sort(unique(stim(istim,:)));
        for iuni = 1:length(uniq);
            un_imP(:,:,iuni) = mean(imP(:,:,stim(istim,:)==uniq(iuni)),3);
        end        
        if length(uniq)>2
            un_imP2 = [];
            un_imP2(:,:,1) = mean(un_imP(:,:,1:floor(length(uniq)/2)),3);
            un_imP2(:,:,2) = mean(un_imP(:,:,floor(length(uniq)/2)+1:end),3);
            un_imP = un_imP2;        
        end
        
        [~,imx] = max(un_imP,[],3);
        
        map(:,:,istim) =rot90(imx',2);
        
%         subplot(2,2,istim)
%         imagesc(map(:,:,istim))
%         if istim==1
%             title('Azimuth')
%                subplot(2,2,istim+2)
%         imagesc(uniq)
%         else
%             title('Elevation')
%                subplot(2,2,istim+2)
%         imagesc(uniq')
%         end
%      
%         title('Monitor')
    end

    a = map(:,:,1)==1&map(:,:,2)==1;
    b = map(:,:,1)==1&map(:,:,2)==2;
    c = map(:,:,1)==2&map(:,:,2)==1;
    d = map(:,:,1)==2&map(:,:,2)==2;
    
    r = cat(3,1,0,0);
    y = cat(3,1,1,0);
    bb = 1-y;
    g = cat(3,0,1,0);
  
    im = bsxfun(@times,a,g)+bsxfun(@times,b,bb)+bsxfun(@times,c,r)+bsxfun(@times,d,y);
    figure
    subplot(121)
    imshow(imrotate(im,-30))
    subplot(122)
    im = imread(getLocalPath('/lab/users/Manolis/Matlab/Datajoint2/+vis2p/@OptImage/map.png'));
    imshow(im)
%     m = cat(3,1,0,1);
%     y = cat(3,1,1,0);
%     b = 1-y;
%     g = 1-m;
%     amp = imP;
%     v = imA;
%     img = bsxfun(@times, amp(:,:,1), b) ...
%         + bsxfun(@times, amp(:,:,2), m) ...
%         + bsxfun(@times, amp(:,:,3), y) ...
%         + bsxfun(@times, size(:,:,4), g);
%     
%     R = img(:,:,1);
%     G = img(:,:,2);
%     B = img(:,:,3);
%     R = R - median(R(:));    R = R / quantile(abs(R(:)),0.99);
%     G = G - median(G(:));    G = G / quantile(abs(G(:)),0.99);
%     B = B - median(B(:));    B = B / quantile(abs(B(:)),0.99);
%     
%     img = cat(3,R,G,B)/2+0.5;
%     img = max(0, min(1, img));
%     img= rgb2hsv(img);
    
%     subplot(131)
%     img1 = img;
%     
%     % img1(:,:,2) = normalize(v);
%     img1(:,:,3) = normalize(vessels);
%     img1 = hsv2rgb(img1);
%     imshow(img1)
%     if params.reverse
%         set(gca,'xdir','reverse')
%     end
%     subplot(132)
%     img1 = img;
%     
%     img1(:,:,3) = normalize(v);
%     img1(:,:,2) = normalize(vessels);
%     
%     img1 = hsv2rgb(img1);
%     imshow(img1)
%     if params.reverse
%         set(gca,'xdir','reverse')
%     end
%     subplot(133)
%     img1 = img;
%     
%     img1(:,:,2) = normalize(vessels);
%     img1 = hsv2rgb(img1);
%     imshow(img1)
    if params.reverse
        set(gca,'xdir','reverse')
    end
    if ikey ~= length(keys)
        pause
    end
end