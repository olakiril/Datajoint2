function makeTuples( obj )

setups = fetchn(Experiments,'setup');
setups = unique(setups);

for isetup = 1:length(setups)
    setup = setups{isetup};
    if strcmp(setup,'2P1')
        start = 'Y:/2Pcalibration/at-photon1.neusc.bcm.tmc.edu/';
    else
        start = 'Y:/2Pcalibration/at-s5no1.neusc.bcm.tmc.edu/';
    end
    
    % Check if both AOD and MPScan are used
    scanProgs = dir(start);
    for iscan_prog = 3:length(scanProgs)
        
        %correct for wrongly naming MPScan  as Galvos..
        if strcmp(scanProgs(iscan_prog).name,'AOD')
            tuple.scan_prog = scanProgs(iscan_prog).name;
        else
            tuple.scan_prog = 'MPScan';
        end
        
        % find all different lenses used
        lenses = dir([start scanProgs(iscan_prog).name '/*']);
        for ilens =  3:length(lenses)
            tuple.lens = str2double(lenses(ilens).name);
            
            % find all different wavelengths used
            wvls = dir([start scanProgs(iscan_prog).name '/' lenses(ilens).name '/*']);
            for iwvls = 3:length(wvls)
                tuple.laser_wavelength = str2double(wvls(iwvls).name);
                
                % find all different gdd corrections made
                gdd = dir([start scanProgs(iscan_prog).name '/' lenses(ilens).name '/' ...
                    wvls(iwvls).name '/*']);
                for igdd = 3:length(gdd)
                    tuple.gdd = str2double(gdd(igdd).name);
                    
                    % find all different measurments
                    files = dir([start scanProgs(iscan_prog).name '/' lenses(ilens).name '/' ...
                        wvls(iwvls).name '/'  gdd(igdd).name '/*']);
                    for ifiles = 3:length(files)
                        
                        % get it into sql timestamp form
                        date = files(ifiles).name(1:10);
                        time = files(ifiles).name(12:end-4);
                        time(strfind(time,'-')) = ':';
                        tuple.rec_timestamp = [date ' ' time];
                        
                        % check if key already exists
                        if isempty(LaserPower(tuple))
                            tuple.calibration = importdata([start ...
                                scanProgs(iscan_prog).name '/' lenses(ilens).name '/' ...
                                wvls(iwvls).name '/'  gdd(igdd).name '/' files(ifiles).name]);
                            insert( obj, tuple );
                            tuple = rmfield(tuple,'calibration');
                        end
                    end
                end
            end
        end
    end
end