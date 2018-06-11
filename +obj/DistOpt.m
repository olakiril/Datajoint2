%{
obj.DistOpt (lookup)                                              # Optional parameters for decoding
dist_opt                 : smallint unsigned                         # decoding option index
---
brief="fill out"        : varchar(127)                              # short description, to be displayed in menus
obj_set=""            : varchar(2000)                             # training set groupA-1,groupA-2;groupB groupC;groupD-1,groupD-2
binsize=500             : float                                     # time window in ms to compute the response
repetitions=10          : tinyint                                   # trial grouping
process="yes"           : enum('no','yes')                          # do it or not
select_method="all"     : enum('all','subsample','single','rf')     # cell selection 
units=10                : tinyint                                   # crossvalidation fold
%}

classdef DistOpt < dj.Lookup
end