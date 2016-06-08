%{
acq.AodVolume (manual) # Aod volume recordings

-> acq.Sessions
aod_volume_start_time: bigint           # start session timestamp
---
aod_volume_stop_time=null   : bigint                        # end of session timestamp
aod_volume_filename         : varchar(255)                  # path to the ephys data
x_coordinate                : double                        # X coordinate of the manipulator
y_coordinate                : double                        # Y coordinate of the manipulator
z_coordinate                : double                        # Z coordinate of the manipulator
objective                   : enum('20x')                   # Objective used
x_range                     : double                        # X volume range
y_range                     : double                        # Y volume range
z_range                     : double                        # Z volume range
x_resolution                : double                        # X volume resolution
y_resolution                : double                        # Y volume resolution
z_resolution                : double                        # Z volume resolution
reps                        : int unsigned                  # Z volume resolution
pmt_green                   : double                        # Setting for green PMT
pmt_red                     : double                        # Setting for red PMT
preamplifier_green          : enum('10e-5','10e-6','10e-7') # preamp setting
preamplifier_red            : enum('10e-5','10e-6','10e-7') # preamp setting
attenuator_degrees          : double                        # Attenautor setting
scan_power                  : double                        # Power out of object in mW
depth                       : double                        # Depth below surface in microns
gdd                         : double                        # GDD compensation
%}


classdef AodVolume < dj.Relvar
    properties(Constant)
        table = dj.Table('acq.AodVolume');
    end
    
    methods 
        function self = AodVolume(varargin)
            self.restrict(varargin{:})
        end
        
        function fn = getFileName(self)
            % Return name of data file matching the tuple in relvar self.
            %   fn = getFileName(self)
            aodPath = fetch1(self, 'aod_volume_filename');
            fn = findFile(RawPathMap, aodPath);
        end
        
        function br = getFile(self)
            % Open a reader for the ephys file matching the tuple in self.
            %   br = getFile(self)
            br = aodReader(getFileName(self), 'Volume');
        end

        function time = getHardwareStartTime(self)
            % Get the hardware start time for the tuple in relvar
            %   time = getHardwareStartTime(self)
            cond = sprintf('ABS(timestamper_time - %ld) < 5000', fetch1(self, 'aod_volume_start_time'));
            rel = acq.SessionTimestamps(cond) & acq.TimestampSources('source = "AOD"') & (acq.Sessions * self);
            time = acq.SessionTimestamps.getRealTimes(rel);
        end
        
        function updateVolumeInformation(self)
            av = fetch(self);
            
            for i = 1:length(av)
                tuple = av(i);
                avr = getFile(acq.AodVolume(tuple));
                
                attr.x_range = range(avr.x);
                attr.y_range = range(avr.y);
                attr.z_range = range(avr.z);
                
                attr.x_resolution = length(avr.x);
                attr.y_resolution = length(avr.y);
                attr.z_resolution = length(avr.z);
                
                f = fields(attr);
                for j = 1:length(f)
                    val = fetch1(self & tuple, f{j});
                    if val == 0
                        update(self & tuple, f{j}, attr.(f{j}));
                    elseif val == attr.(f{j})
                    else
                        disp('Value set but does not match file');
                    end
                end
            end
        end

    end
end
