%{
# 
dec_opt                     : smallint unsigned             # decoding option index
---
brief="fill out"            : varchar(127)                  # short description, to be displayed in menus
train_set                   : varchar(2000)                 # training set groupA-1,groupA-2;groupB groupC;groupD-1,groupD-2
test_set                    : varchar(2000)                 # testing set
decoder="fitclinear"        : enum('fitclinear','fitcecoc','fitcsvm') # decoding method
repetitions=10              : tinyint                       # trial grouping
process="yes"               : enum('no','yes')              # do it or not
select_method="all"         : enum('all','subsample','expand') # cell selection
shuffle=3000                : mediumint                     # chance performance
k_fold=10                   : tinyint                       # crossvalidation fold
dec_params=null             : varchar(1024)                 # 
neurons=null                : mediumint                     # 
fold_selection=null         : enum('random','rotation','tilt','scale','light2_ene','y','light1_xloc') # 
%}

classdef DecodeOpt < dj.Lookup
end
