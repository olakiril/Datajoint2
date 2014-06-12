function  [iH iS iV] = plot(obj,varargin)

% plot(obj)
%
% Plots the intrinsic imaging data aquired with Intrinsic Imager
%
% MF 2012-09

params.sigma = 2; %sigma of the gaussian filter
params.exp = 1; % exponent factor of rescaling
params.reverse = 0; % reverse the axis
params.range = 3.14/2; %angle limit

params = getParams(params,varargin);

keys = fetch(obj);

for ikey = 1:length(keys)
      
    [imP vessels imA] = fetch1(vis2p.OptImageBar(keys(ikey)),'ang','vessels','amp');
    
    imA(imA>prctile(imA(:),99)) = prctile(imA(:),99);
    
    [h1 h2] = hist(reshape(imP(imP~=0),[],1),100);
    mxv = h2(h1 == max(h1));
    imP = imP - mxv(1);
    imP(imP<-3.14) = imP(imP<-3.14) +3.14*2;
    imP(imP>3.14) = imP(imP>3.14) -3.14*2;
    imP(imP<0) = -exp((imP(imP<0)+ params.range)*params.exp);
    imP(imP>0) = exp((abs(imP(imP>0)- params.range))*params.exp);
    
    
    h = normalize(imP);
    s = ones(size(imP));
    v = normalize(imA);
    s2 = normalize(imA);
    v2 = normalize(vessels);
    
    if nargout>0
        iH{ikey} = h;
        iS{ikey} = s2;
        iV{ikey} = v2;
    else
        figure
        set(gcf,'position',[50 200 920 435])
        set(gcf,'name',['OptMap ' keys(ikey).exp_date ' ' num2str(keys(ikey).scan_idx)])
        
        subplot(121)
        im = (hsv2rgb(cat(3,h,cat(3,s,v))));
        im = imgaussian(im,params.sigma);
        imshow(im)
        if params.reverse
            set(gca,'xdir','reverse')
        end
        
        subplot(122)
        im = (hsv2rgb(cat(3,h,cat(3,s2,v2))));
        im = imgaussian(im,params.sigma);
        imshow(im)
        
        if params.reverse
            set(gca,'xdir','reverse')
        end
        if ikey ~= length(keys)
            pause
        end   
    end
end

if ikey == 1 & nargout>0
    iH = iH{1};
    iS = iS{1};
    iV = iV{1};
end