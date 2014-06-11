% Scans/play
% function play( obj, frames )
%
%  plays the two-photon movie in real time or faster 
%  The optional parameter 'frames' lists the frames to play


function play( obj, frames )

for key = enumerate( obj )
    disp('Playing movie for');
    disp(key);
    tpr = tpReader(Scans(key));
    gChan = getImagingChannel( tpr, 1 );
    rChan = getImagingChannel( tpr, 2 );
    sz = size(gChan);

    if nargin<=1
        frames=1:sz(3);
    end

    % play
    [rasterCorrection,motionCorrection] = take( Movies(key), 'raster_correction', 'motion_correction' );
    for iFrame = 1:length(frames)
        subplot(131);
        img = anscombe( cat( 3, rChan(:,:,frames(iFrame))/2, gChan(:,:,frames(iFrame)), zeros(sz(1:2)) ) );
        imshow( img );

        %M(iFrame) = getframe;
        subplot(132);
        img(:,:,1) = fixInterlace( img(:,:,1), 'initialSolution', rasterCorrection(iFrame,:), 'maxIter', 0 );
        img(:,:,2) = fixInterlace( img(:,:,2), 'initialSolution', rasterCorrection(iFrame,:), 'maxIter', 0 );
        img(:,:,3) = fixInterlace( img(:,:,3), 'initialSolution', rasterCorrection(iFrame,:), 'maxIter', 0 );
        imshow( img );
        
        subplot(133);
        img(:,:,1) = alignFrame( img(:,:,1), motionCorrection(iFrame,:));
        img(:,:,2) = alignFrame( img(:,:,2), motionCorrection(iFrame,:));
        img(:,:,3) = alignFrame( img(:,:,3), motionCorrection(iFrame,:));
        imshow( img );
        drawnow;
    end

    %movie(M); 
end




function img = anscombe( img )
   img = sqrt(max(0,double(img)+4))/sqrt(2052);