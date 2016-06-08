%{
olf.TracesOpt (lookup) # 
trace_opt       : smallint unsigned      # 
---
brief="fill out"            : varchar(127)                  # m) short description, to be displayed in menus
trace_computation="raw"     : enum('oopsiLnorm','dfof','calcDFoF','spikeExtraction','rnorm,calcDFoF','rnorm,dfof','gnorm,calcDFoF','gnorm,dfof','gnotm,spikeExtraction','rnorm,spikeExtraction','oopsi','oopsiThr','rnorm,oopsi','rnormMinus,oopsi','rnormMinus,calcDFoF','oopsiM','rtrace','rtrace,oopsi','gnorm,oopsi','antrace','gnormMinus,calcDFoF','gnormMinus,oopsi','antrace,calcDFoF','antrace,oopsi','rnorm,gnorm,calcDFoF','gnormAll,oopsi','pcMinus,calcDFoF','pcMinus,oopsi','pcMinus','oopsiBest','pcMinus,oopsiAll','oopsiAll','pcMinusHighPass,oopsi','pcMinusHighPass,oopsiAll','pcMinusHighPass','gnormMinus','gpcMinus,oopsi','globalPCminus','globalPCminus,oopsi','raw') # c) method of trace computation
highPass=null               : float                         # m) (Hz) stop frequency 1 of highpass filter for baseline computation
lowPass=null                : float                         # m) (Hz) stop frequency 2 of low-pass filter
lamda=null                  : float                         # Average firing rate
tau=null                    : float                         # m) time constant, e.g. for the exponential model, for example.
sigma=null                  : float                         # 
process="yes"               : enum('no','yes')              # m) do it or not
%}


classdef TracesOpt < dj.Relvar
	methods

		function self = TracesOpt(varargin)
			self.restrict(varargin{:})
		end
	end

end