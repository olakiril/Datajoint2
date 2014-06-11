%{
vis2p.StatsOriParams (lookup) # 
ori_param       : tinyint                # m) the parameter settings
---
ori_num=16                  : smallint                      # m) number of different orientations
process="yes"               : enum('no','yes')              # m) compute or not compute
lambda=80                   : smallint                      # m) lambda gabor parameter
psi=0.00                    : float(4,2)                    # m) psi gabor parameter
gamma=0.80                  : float(4,2)                    # m) gamma gabor parameter
bw=1.00                     : float(4,2) unsigned zerofill  # m) bw gabor parameter
rf_snr_thr=1.0              : float(2,1) unsigned           # m) rf snr threshold
rf_p_thr=0.0500             : float(5,4) unsigned           # m) rf p threshold
%}


classdef StatsOriParams < dj.Relvar
	methods

		function self = StatsOriParams(varargin)
			self.restrict(varargin{:})
		end
	end

end