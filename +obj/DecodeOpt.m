%{
obj.DecodeOpt (lookup)                                              # Optional parameters for decoding
dec_opt                 : smallint unsigned                         # decoding option index
---
brief="fill out"        : varchar(127)                              # short description, to be displayed in menus
train_set=""            : varchar(2000)                             # training set groupA-1,groupA-2;groupB groupC;groupD-1,groupD-2
test_set=""             : varchar(2000)                             # testing set
binsize=500             : float                                     # time window in ms to compute the response
decoder="fitclinear"    : enum('fitclinear','fitcecoc','fitcsvm')   # decoding method
repetitions=10          : tinyint                                   # trial grouping
trial_method="random"   : enum('random','sequential')               # trial selection method
process="yes"           : enum('no','yes')                          # do it or not
rf_opt=0                : smallint unsigned                         # restrict anaysis to population RF
select_method="all"     : enum('all','subsample','single','rf')     # cell selection 
shuffle=1000            : mediumint                                 # chance performance
k_fold=10               : tinyint                                   # crossvalidation fold
%}

classdef DecodeOpt < dj.Relvar
	methods
		function self = DecodeOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end