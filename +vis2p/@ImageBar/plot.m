function  plot(obj,varargin)

% plot(obj)
%
% Plots the intrinsic imaging data aquired with Intrinsic Imager
%
% MF 2012-09

params.sigma = 2; %sigma of the gaussian filter
params.exp = 1; % exponent factor of rescaling
params.reverse = 0; % reverse the axis

params = getParams(params,varargin);

figure
set(gcf,'position',[50 200 920 435])

keys = fetch(obj);

for ikey = 1:length(keys)
    
    set(gcf,'name',['OptMap ' keys(ikey).exp_date ' ' num2str(keys(ikey).scan_idx)])
    
    [imP vessels imA] = fetch1(ImageBar(keys(ikey)),'ang','vessels','amp');
    
    imA(imA>prctile(imA(:),99)) = prctile(imA(:),99);
    
    
    [h1 h2] = hist(reshape(imP(imP~=0),[],1),100);
    mxv = h2(h1 == max(h1));
    imP = imP - mxv;
    imP(imP<-3.14) = imP(imP<-3.14) +3.14*2;
    imP(imP>3.14) = imP(imP>3.14) -3.14*2;
    s = 3;
    imP(imP<0) = -exp((imP(imP<0)+ 3.14/2)*params.exp);
    imP(imP>0) = exp((abs(imP(imP>0)- 3.14/2))*params.exp);
    
    subplot(121)
    h = normalize(imP);
    s = ones(size(imP));
    v = normalize(imA);
    im = (hsv2rgb(cat(3,h,cat(3,s,v))));
    im = imgaussian(im,params.sigma);
    imshow(im)
    if params.reverse
        set(gca,'xdir','reverse')
    end
    
    subplot(122)
    s = normalize(imA);
    v = normalize(vessels);
    im = (hsv2rgb(cat(3,h,cat(3,s,v))));
    im = imgaussian(im,params.sigma);
    imshow(im)
    
    if params.reverse
        set(gca,'xdir','reverse')
    end
    
    if ikey ~= length(keys)
        pause
    end
end