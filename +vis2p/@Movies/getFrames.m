function frames = getFrames(self, channel, frameIdx, raster, mc)

% function frames = getFrames(self, channel, frameIdx, raster, warp, wDeg)
%
% gets the frames of a scan and applies corrections if provided
%
% MF 2013-08

if isobject(self)
    % get scan
    scim = tpReader(Scans.*self);
else
    scim = self;
end

% chech for input parameters
if nargin<5; mc = fetch1(self, 'motion_correction');
    motion = mc;
if nargin<4; raster = fetch1(self,'raster_correction');
    if isnan(raster); raster = false; end
if nargin<3; frameIdx = 1:fetch1(self,'nframes');
if nargin<2; channel = 1;end;end;end;end %#ok<*ALIGN>

%  motion = mc.warp_degree;
%     degrees = [mc.warp_degree mc.warp_degree size(mc.warp_polynom,2)-2*(mc.warp_polynom+1)];
%     if isnan(mc.warp_degree);  mc.warp_degree = false; end;end;end;end;end

% read specified frames
frames = single(scim.read(channel, frameIdx));

% apply raster correction
if raster
    disp 'raster correction...'
    frames = tpMethods.RasterCorrection.apply(frames, raster);
end

% apply motion compansation
if mc
    disp 'motion correction...'
    frames = tpMethods.MotionCorrection.apply(frames, motion);
end
% % apply motion compansation
% if mc.warp_degree
%     disp 'motion correction...'
%     for i=1:length(frameIdx)
%         frames(:,:,i) = tpMethods.YWarp.apply(frames(:,:,i), motion(frameIdx(i),:), degrees);
%     end
% end

