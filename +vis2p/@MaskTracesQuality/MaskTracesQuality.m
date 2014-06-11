%{
vis2p.MaskTracesQuality (computed) # 
-> vis2p.MaskGroup
masknum         : mediumint              # 
---
ca_event_snr                : float                         # c) the ratio of the fitted model stddev to the stddev of the residual
cell_ts=CURRENT_TIMESTAMP   : timestamp                     # c) automatic timestamp
%}


classdef MaskTracesQuality < dj.Relvar
	methods

		makeTuples(self, key)
		

		function self = MaskTracesQuality(varargin)
			self.restrict(varargin{:})
		end
	end

end