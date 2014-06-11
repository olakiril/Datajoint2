%{
vis2p.StatsImagesParams (lookup) # 
statparam       : smallint unsigned      # 
---
imresize=1                  : float                         # m) image resize ratio
wrf="yes"                   : enum('no','yes')              # m) only include pixels within the RF
matrix="no"                 : enum('[8 10]','[2 3]','no')   # 
process="yes"               : enum('no','yes')              # m) do it or not
%}


classdef StatsImagesParams < dj.Relvar
	methods

		function self = StatsImagesParams(varargin)
			self.restrict(varargin{:})
		end
	end

end