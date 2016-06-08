%{
acq.AodScan (manual) # Aod scans

-> acq.Sessions
aod_scan_start_time: bigint             # start session timestamp
---
aod_scan_stop_time=null     : bigint                        # end of session timestamp
aod_scan_filename           : varchar(255)                  # path to the ephys data
x_coordinate                : double                        # X coordinate
y_coordinate                : double                        # Y coordinate
z_coordinate                : double                        # Z coordinate
objective                   : enum('20x')                   # Objective used
motion_planes               : tinyint                       # 
pmt_green                   : double                        # Setting for green PMT
pmt_red                     : double                        # Setting for red PMT
preamplifier_green          : enum('10e-5','10e-6','10e-7') # preamp setting
preamplifier_red            : enum('10e-5','10e-6','10e-7') # preamp setting
attenuator_degrees          : double                        # Attenautor setting
scan_power                  : double                        # Power out of object in mW
depth                       : double                        # Depth below surface in microns
gdd                         : double                        # GDD compensation
%}


classdef AodScan < dj.Relvar
    properties(Constant)
        table = dj.Table('acq.AodScan');
    end
    
    properties(Constant,Access=public)
        x_step = 100 / 146000000;
        y_step = 100 / 146000000;
        z_step = 100 / 70000;
    end
    
    methods 
        function self = AodScan(varargin)
            self.restrict(varargin{:})
        end
        
        function fn = getFileName(self)
            % Return name of data file matching the tuple in relvar self.
            %   fn = getFileName(self)
            aodPath = fetch1(self, 'aod_scan_filename');
            fn = findFile(RawPathMap, aodPath);
        end
        
        function br = getFile(self, dataType)
            % Open a reader for the ephys file matching the tuple in self.
            %   br = getFile(self)
            
            if nargin < 2
                dataType = 'Functional';
            end
            
            br = aodReader(getFileName(self), dataType);
        end
        
        function time = getHardwareStartTime(self)
            % Get the hardware start time for the tuple in relvar
            %   time = getHardwareStartTime(self)
            cond = sprintf('ABS(timestamper_time - %ld) < 5000', fetch1(self, 'aod_scan_start_time'));
            rel = acq.SessionTimestamps(cond) & acq.TimestampSources('source = "AOD"') & (acq.Sessions * self);
            time = acq.SessionTimestamps.getRealTimes(rel);
        end
    end
end
