function makeTuples( obj, key )

tuple = key;

% Load image
field = fetch1(vis2p.RFMap(key),'onpoff_rf');

% fit gauss
s = fetch1(vis2p.VisStims(key),'stim_file');
theta = atand(s.params.constants.monitorSize(1)/2/s.params.constants.monitorDistance)*2;
deg2dot = s.params.constants.dotNumX/theta;
tuple.gauss_fit = fitGauss(field,'deg2dot',deg2dot);

% compute snr
tuple.snr = computeRFsnr(tuple.gauss_fit,field);

% correct for nan values
tuple.snr(isnan(tuple.snr) | isinf(tuple.snr)) = 0;

% insert things
insert( obj, tuple );




