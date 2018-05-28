%{
# 
rdm_opt       : smallint unsigned      # 
---
->obj.ReDimFunc
brief="fill out"            : varchar(127)                        # short description, to be displayed in menus
binsize=500                 : float                               # time window in ms to compute the response
dimensions=3                : mediumint                           # final dimensions
gauss_win=1500              : int                                 # trial time bluring in msec
select_method="all"         : enum('all','subsample','rf')           # cell selection 
cell_num=50                 : int                                 # cell number for restricted selections 
process="yes"               : enum('no','yes')                    # do it or not
%}


classdef ReDimOpt < dj.Lookup
end