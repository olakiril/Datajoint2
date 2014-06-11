%{
vis2p.MaskCells (imported) # 
mouse_id        : smallint unsigned      # 
exp_date        : date                   # 
scan_idx        : smallint unsigned      # 
masknum         : mediumint              # 
---
img_x                       : decimal(11,4)                 # i) pixel coordinate of the cell's center
img_y                       : decimal(11,4)                 # i) pixel's coordinate of the cell's center
img_z=0.0000                : decimal(11,4)                 # i) pixel's coordinate of the cell's center 
cell_radius                 : decimal(6,4) unsigned         # i) cell radius in pixels 
green_contrast              : decimal(6,4)                  # i) cell contrast (natural log scale)
red_contrast=null           : decimal(6,4)                  # i) cell contrast in the red image (gamma scale)
sharpness                   : decimal(5,2) unsigned         # i) cell sharpness
fit_qual                    : decimal(4,3) unsigned         # i) correlation between fitted model and data.  Should be > 0.95
%}


classdef MaskCells < dj.Relvar
	methods

		makeTuples(self, key)
	

		function self = MaskCells(varargin)
			self.restrict(varargin{:})
		end
	end

end