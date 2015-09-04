%{
vis2p.Movies (imported) #
-> vis2p.Scans
---
xsize                       : smallint unsigned             # i) the x frame size of the scan
ysize                       : smallint unsigned             # i) the y frame size of the scan
zsize=1                     : smallint unsigned             # i) the z frame size of the scan
nframes                     : smallint unsigned             # i) the number of frames
fps                         : double                        # i) frames per second
ephys_fs=null               : double                        # electrophysiology rate
raw_green                   : mediumblob                    # i) mean green image before any corrections
raw_red=null                : mediumblob                    # i) mean red image before any corrections
affine_green=null           : mediumblob                    # i) mean green image after affine correction. No motion correction
affine_red=null             : mediumblob                    # i) mean red image after affine correction. No motion correction
alignXshifts=null           : mediumblob                    # i) the aligning X shift of each frame
alignYshifts=null           : mediumblob                    # i) the aligning Y shift of each frame
alignZshifts=null           : mediumblob                    # i) the aligning Z shift of each frame
motion_correction=null      : mediumblob                    #
raster_correction=null      : mediumblob                    #
%}


classdef Movies < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = vis2p.Scans('aim <> "stack" and problem_type = "none!"') & vis2p.Experiments('process = "yes"')
    end
    
    methods(Access=protected)
                
        makeTuples( obj, key )
    end
    
    methods
        
        function self = Movies(varargin)
            self.restrict(varargin{:})
        end
                
        obj = compute( obj, fieldname, key )
        
        do( obj )
        
        frames = getFrames(obj, channel, frameIdx, raster, mc,scim)

        volume = getAODVolume(obj,key)
        
        plotEye(obj,clim)
        
        plotMovie(obj,tracetype)
        
        plotSite(obj,tracetype)

        tpr = tpReader( obj )
        
    end
    
end