%{
olf.DecodeOpt (lookup)                                              # Optional parameters for decoding
dec_opt                 : smallint unsigned                         # decoding option index
---
discription="fill out"  : varchar(127)                              # short description, to be displayed in menus
train_set=""            : varchar(2000)                             # training set groupA-1,groupA-2;groupB groupC;groupD-1,groupD-2
test_set=""             : varchar(2000)                             # testing set
decoder="fitclinear"    : enum('fitclinear','fitcecoc','fitcsvm')   # decoding method
repetitions=1           : tinyint                                   # trial grouping
trial_method="random"   : enum('random','sequential')               # trial selection method
process="yes"           : enum('no','yes')                          # do it or not
shuffle=1000            : mediumint                                 # chance performance
k_fold=1                : tinyint                                   # crossvalidation fold
%}

classdef DecodeOpt < dj.Relvar
	methods
		function self = DecodeOpt(varargin)
			self.restrict(varargin{:})
		end
    end
end