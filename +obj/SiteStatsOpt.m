%{
stats_opt               : smallint unsigned                         # stats option index
---
brief="fill out"        : varchar(127)                              # short description, to be displayed in menus
binsize=500             : float                                     # time window in ms to compute the response
neurons=128             : smallint                                  # number of neurons to evaluate pop stats
reps=128                : smallint                                  # repetitions for neuron sampling
process="yes"           : enum('no','yes')                          # do it or not
%}

classdef SiteStatsOpt < dj.Lookup
end