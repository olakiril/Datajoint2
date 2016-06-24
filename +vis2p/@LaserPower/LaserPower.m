%{
vis2p.LaserPower (imported) # 
rec_timestamp   : timestamp              # i) the timestamp..
scan_prog       : enum('MPScan','AOD')   # i) Scanning program
lens            : tinyint unsigned       # i) lens magnification
laser_wavelength: smallint unsigned      # i) (nm)
gdd             : mediumint unsigned     # i) gdd correction
---
calibration                 : mediumblob                    # i) the measurments.
%}


classdef LaserPower < dj.Relvar 

% 	properties
% 		popRel = vis2p.LaserPower
% 	end

	methods(Access=protected)
		makeTuples( obj )
	end

end