function [dataStruct] = LoadAFMData(pathname,varargin)
%LOADAFMDATA Load AFM Data from an External Source
%   This function takes in a pathname, and will then parse it looking for
%   relevant AFM SFS files. Implementing new filetypes involves adding new
%   cases to the switch-case that ties the label (i.e. first string in the
%   filename, which is underscore "_" separated) to the correct steps for
%   processing. For example, if you are attempting to load
%   AFMData_Experiment1.ibw, then there must be a case for when strSwitch
%   is "AFMData". That case should (for the first switch-case) determine
%   which files in the directory are the correct ones to load. Then, the
%   second switch case requires the user to write the code required to load
%   at least the following variables: z (Z-Sensor [m]) & d (Deflection [m]).
%   
%   The file will then organize and clean up the input data, find the
%   region of interest for fitting, and also perform averaging operations
%   for all of the files that roughly share an approach velocity. The
%   original data is kept in the first 1-N rows (where N is the number of
%   files provided), and the new "averaged" datasets are stored at the END
%   of dataStruct. There is one row for each approach velocity detected. 

% Default values
includeRetract = 0;             % Include data from the retract curve
filterType = 'none';            % Choose the filter used to smooth data
N = 2;                          % Order of Butterworth filter, if used
cutoff_Hz = 5000;               % Cutoff frequency
findRep = 'forward';            % Search direction for the repulsive region
removeNegatives = 1;            % Remove negative values in the data stream
useBruker = 0;                  % Use the Bruker Toolbox to load .spm files

% Read varargin values
if ~isempty(varargin)
    if ~isa(varargin{1},'struct')
        error('You can only pass a STRUCT containing settings to LoadAFMData(). Please verify you have done so, because the type does not appear to be "struct".');
    else
        % Load the settings
        inputSettings = varargin{1};
        
        % Load mandatory settings
        includeRetract = inputSettings.includeRetract;
        filterType = inputSettings.filterType;
        findRep = inputSettings.findRep;
        removeNegatives = inputSettings.removeNegatives;
        useBruker = inputSettings.useBruker;
        
        % Load conditional settings
        if strcmpi(filterType,'butter')
            N = inputSettings.N;
            cutoff_Hz = inputSettings.cutoff_Hz;
        end
    end
end

FilesCheck=dir([pathname '/*.*']);

% Remove Directories
FilesCheck=FilesCheck(~ismember({FilesCheck.name},{'.','..'}));
toRemove = find([FilesCheck.isdir] == 1);
FilesCheck(toRemove) = [];

% Remove Filetypes To Ignore
toRemove = find(~endsWith({FilesCheck.name}, {'.ibw','.txt','.spm','.mat','.csv'}));
FilesCheck(toRemove) = [];

toRemove = find(contains({FilesCheck.name}, {'settingsStruct','Settings'}));
FilesCheck(toRemove) = [];

for i = 1:length(FilesCheck)
    FilesTempPrep = FilesCheck(i).name;
    FilesTempPrep = strsplit(FilesTempPrep, {'_' '.'},'CollapseDelimiters',true);
    strSwitchPrep = FilesTempPrep{1};
    tempLabels{i} = strSwitchPrep;
end

if exist('tempLabels','var')
    if length(unique(tempLabels)) > 1
        error('ERROR: Attempted to pass previous fit results as input. Please ensure only data is stored in this directory.');
    end
else
    error('ERROR: No valid AFM files found! Please ensure you have selected the correct directory, and that LoadAFMData contains a case describing how to load that filetype.');
end

if length(FilesCheck) > 1
    FilesTemp = FilesCheck(1).name;
    FilesTemp = strsplit(FilesTemp, {'_' '.'},'CollapseDelimiters',true);
    
    strSwitch = FilesTemp{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            Files=dir([pathname '/*.ibw']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            end
        
        case lower({'CytoDCell-MelanocytesExp', 'UCell-MelanocytesExp',...
                'CytoDcell-MelanomaExp', '1205LuWT-MelanomaExp',...
                'U1-MelanomaExp', 'U2-MelanomaExp', 'HFFD1-FibroblastsExp',...
                'HFFD2-FibroblastsExp'})
            Files=dir([pathname '/*.spm']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
                
                experiment_number(k) = str2num(strrep(FileInfo{2},'-','.'));    % Experiment Number
                point_number(k) = str2num(strrep(FileInfo{3},'-','.'));         % Cell Number
                run_number(k) = str2num(strrep(FileInfo{4},'-','.'))+1;         % Run Number
            end
        
    end
    
else
    
    FilesTemp = FilesCheck.name;
    FilesTemp = strsplit(FilesTemp, {'_' '.'},'CollapseDelimiters',true);
    
    strSwitch = FilesTemp{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            Files=dir([pathname '/*.ibw']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number

        case lower({'CytoDCell-MelanocytesExp', 'UCell-MelanocytesExp',...
                'CytoDcell-MelanomaExp', '1205LuWT-MelanomaExp',...
                'U1-MelanomaExp', 'U2-MelanomaExp', 'HFFD1-FibroblastsExp',...
                'HFFD2-FibroblastsExp'})
            Files=dir([pathname '/*.spm']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
                
            experiment_number = str2num(strrep(FileInfo{2},'-','.'));    % Experiment Number
            point_number = str2num(strrep(FileInfo{3},'-','.'));         % Cell Number
            run_number = str2num(strrep(FileInfo{4},'-','.'))+1;   
           
    end
end

% Initialize our output structure
dataStruct = struct;

% Remove the files we don't care about
FilesRemove=(~ismember({Files.name},{FilesCheck.name}));
Files(FilesRemove) = [];

% Loop through the files and proceed with loading, parsing, and calculation
% for each one.
for k = 1:length(Files)
    FileNames = Files(k).name;
    FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
    
    strSwitch = FileInfo{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            % This case shows the steps used to load .ibw files from an
            % MFP3D running the Asylum AFM software. The file header is
            % read into headerValue and all of the relevant settings are
            % extracted from there. The actual data itself relies on the
            % function IBWread to load the AFM experiment observables. Note
            % that IBWread() was NOT created by the authors of this
            % repository/manuscript. A license for that script is included
            % in this directory, and named IBWread_license.txt. Credit goes
            % to Jakub Bialek (2009).
            RawData = IBWread([Files(k).folder '/' Files(k).name]);
            [headerValue,~] = strsplit(RawData.WaveNotes,'\r',...
            'DelimiterType','RegularExpression');

            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, m/s
            point_number_temp = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number_temp = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            
            dataStruct(k).z = RawData.y(:,1);
            dataStruct(k).d = RawData.y(:,2);
                        
            % Find Spring Constant
            varIndex = find(contains(headerValue,'SpringConstant'),1);
            temp = headerValue{varIndex};
            temp(isspace(temp)) = [];
            temp = split(temp,':');
            k_cantilever(k) = str2num(temp{2}); % Value in N/m
            dataStruct(k).k_cantilever = k_cantilever(k);
            
            % Find Deflection InvOLS
            varIndex = find(contains(headerValue,'InvOLS'));
            temp = headerValue{varIndex};
            temp(isspace(temp)) = [];
            temp = split(temp,':');
            defl_InVOLS(k) = str2num(temp{2}); % Value in nm/V
            
            % Create time array
            varIndex = find(contains(headerValue,'NumPtsPerSec'));
            temp = headerValue{varIndex};
            temp(isspace(temp)) = [];
            temp = split(temp,':');
            dataStruct(k).dt = 1/str2num(temp{2}); % Value in s
            dataStruct(k).n_sam = RawData.Nsam;
            dataStruct(k).t = dataStruct(k).dt.*((1:size(dataStruct(k).d))-1)';
    
        case lower({'CytoDcell', '1205LuWT', 'CytoDCell-MelanocytesExp',...
                'UCell-MelanocytesExp', 'CytoDcell-MelanomaExp',...
                '1205LuWT-MelanomaExp', 'U1-MelanomaExp', 'U2-MelanomaExp',...
                'HFFD1-FibroblastsExp', 'HFFD2-FibroblastsExp'})
            
            if useBruker
                NSMU = NSMatlabUtilities();

                %open a file
                NSMU.Open(fullfile([Files(k).folder '/' Files(k).name]));

                %Get Deflection Sensitivity and Units
                if NSMU.IsOpen
                    defl_InVOLS = NSMU.GetDeflSensitivity();
                    deflSensUnits = NSMU.GetDeflSensitivityUnits();
                    if (bDisplay)
                        disp ' ';
                        disp(['Deflection Sensitivity: ', num2str(defl_InVOLS), ' ' deflSensUnits]);
                    end
                end
                
                %Get F vs Z plot of channel 1
                [zApproach, zRetract, FApproach, FRetract, xLabel, yLabel] = NSMU.CreateForceZPlot(1, NSMU.FORCE, 0);

                zString = string(regexp(xLabel,'(?<=\()\S+(?=\))','match'));
                switch zString
                    case 'nm'
                        zScale = 1e-3;
                    case 'um'
                        zScale = 1;
                    otherwise
                        zScale = 1;
                end

                zApproach = zApproach.*zScale;
                zRetract = zRetract.*zScale;

                point_number_temp = str2num(strrep(FileInfo{3},'-','.'));     % Cell/Position Number
                run_number_temp = str2num(strrep(FileInfo{4},'-','.'))+1;     % Run Number
                v_approach(k) = NSMU.GetForwardRampVelocity(1, NSMU.METRIC).*1e-9;  % Approach Velocity, m/s

                dataStruct(k).t = time;
                dataStruct(k).dt = mean(diff(time));
                dataStruct(k).n_sam = length(time);

                k_cantilever(k) = NSMU.GetSpringConstant();
                dApproach = (FApproach./k_cantilever(k));
                dRetract = (FRetract./k_cantilever(k));

                dataStruct(k).z = [zApproach;zRetract]*(1e-6);
                dataStruct(k).d = [dApproach;dRetract]*(1e-9);

                NSMU.Close();
                
            else
               
                fid = fopen([Files(k).folder '/' Files(k).name],'r');

                sHSDC_datalength=di_header_find([Files(k).folder '/' Files(k).name],'\Data length:');
                sHSDC_dataoffset=di_header_find([Files(k).folder '/' Files(k).name],'\Data offset:');

                for i_loop = 1:length(sHSDC_datalength)
                    fseek(fid,sHSDC_datalength(i_loop),-1);
                    line = fgets(fid);
                    datalength(i_loop) = extract_num(line);

                    if i_loop > 1
                        fseek(fid,sHSDC_dataoffset(i_loop-1),-1);
                        line = fgets(fid);
                        dataoffset(i_loop-1) = extract_num(line);
                    end
                end

                frewind(fid);
                fseek(fid,datalength(1),-1);
                rawData = (fread(fid, datalength(1)/2, 'int16'));

                usefulData = rawData;

                a_ind = [1 datalength(2)/4+1 datalength(2)/2+1 datalength(2)/2+datalength(3)/4+1];
                b_ind = [datalength(2)/4 datalength(2)/2 (datalength(2)/2)+(datalength(3)/4) (datalength(2)/2+datalength(3)/2)];
                
                if b_ind(end) > length(usefulData)
                    b_ind(end) = length(usefulData);
                end
                
                approachcurvesraw=flip(usefulData(a_ind(1):b_ind(1)));
                retractcurvesraw=flip(usefulData(a_ind(2):b_ind(2)));
                zapproachraw=flip(usefulData(a_ind(3):b_ind(3)));
                zretractraw=flip(usefulData(a_ind(4):b_ind(4)));

                k_temp = di_header_find([Files(k).folder '/' Files(k).name],'\Spring Constant:');
                fseek(fid,k_temp(1),-1);
                line = fgets(fid);
                k_cantilever(k) = extract_num(line);

                InVOLS_temp = di_header_find([Files(k).folder '/' Files(k).name],'\@Sens. DeflSens:');
                fseek(fid,InVOLS_temp(1),-1);
                line = fgets(fid);
                defl_InVOLS = extract_num(line); % Value in nm/V

                InVOLS2_temp = di_header_find([Files(k).folder '/' Files(k).name],'\@Sens. ZsensSens:');
                fseek(fid,InVOLS2_temp(1),-1);
                line = fgets(fid);
                zSens_InVOLS = extract_num(line); % Value in nm/V

                deflScale_temp = di_header_find([Files(k).folder '/' Files(k).name],'\@4:Z scale: V [Sens. DeflSens]');
                fseek(fid,deflScale_temp,-1);
                line = fgets(fid);
                SensDeflSens_InVOLS = extract_num(line); % Value in V/LSB

                zScale_temp = di_header_find([Files(k).folder '/' Files(k).name],'\@4:Z scale: V [Sens. ZsensSens]');
                fseek(fid,zScale_temp,-1);
                line = fgets(fid);
                SensZSens_InVOLS = extract_num(line); % Value in V/LSB

                deflScaling = defl_InVOLS*SensDeflSens_InVOLS;
                zPosnScaling = zSens_InVOLS*SensZSens_InVOLS;

                approachDefl = (approachcurvesraw(:)-(approachcurvesraw(1))).*deflScaling;
                retractDefl = (retractcurvesraw(:)-(approachcurvesraw(1))).*deflScaling;
                zApproach = (zapproachraw(:)-(zapproachraw(1))).*zPosnScaling.*(1e6).*(1e-9);
                zRetract = (zretractraw(:)-(zapproachraw(1))).*zPosnScaling.*(1e6).*(1e-9);
                
                % Make time array
                frewind(fid);
                dt_temp = di_header_find([Files(k).folder '/' Files(k).name],'\@Sample period: V (0.1000000 us/LSB)');
                fseek(fid,dt_temp(2),-1);
                line = fgets(fid);
                line = erase(line,'\@Sample period: V (0.1000000 us/LSB)');
                dt = extract_num(line).*1e-6; % Value in us
                time = ((1:size(approachDefl,1)).*dt)';

                point_number_temp = str2num(strrep(FileInfo{3},'-','.'));   % Cell/Position Number
                run_number_temp = str2num(strrep(FileInfo{4},'-','.'))+1;   % Run Number
                
                buffer = 10;
                dataError_pos = ( length(zApproach) - find(flip(gradient(zApproach(1:end-buffer)))<0,1) );    % Find last position above 0
                dataError_pos2 = ( length(approachDefl) - find(abs(flip(gradient(approachDefl)))>10*max(abs(gradient(approachDefl(end-buffer:end)))),1) );    % Find where deflection blows up
                dataError_pos = min([dataError_pos dataError_pos2]);    % Choose worst cutoff

                if isempty(dataError_pos) || dataError_pos == 0
                    dataStruct(k).z = zApproach*(1e-6);
                    dataStruct(k).d = approachDefl*(1e-9);
                    dataStruct(k).t = time;
                    dataStruct(k).dt = dt;
                    dataStruct(k).n_sam = length(time);
                else
                    dataStruct(k).z = zApproach(dataError_pos+buffer:end)*(1e-6);
                    dataStruct(k).d = approachDefl(dataError_pos+buffer:end)*(1e-9);
%                     dataStruct(k).t = time(1:(dataError_pos+buffer)+1);
                    dataStruct(k).t = time(1:length(dataStruct(k).d));
                    dataStruct(k).dt = dt;
                    dataStruct(k).n_sam = length(dataStruct(k).t);
                end
                
                v_approach(k) = mean(abs(diff(dataStruct(k).z))./(dataStruct(k).dt*10));    % Approach Velocity, m/s
                fclose(fid);
                
            end
            
    end
    
    % Save velocity
    dataStruct(k).v_approach = v_approach(k);
    dataStruct(k).point_number = point_number(k);
    dataStruct(k).run_number = run_number(k);
    
    % Pre-Processing
    % Find the approach portion of the data
    [~, z_max_ind] = max(dataStruct(k).z);
    
    % Make sure we are handling row vectors
    if isrow(dataStruct(k).z) dataStruct(k).z = dataStruct(k).z'; end
    if isrow(dataStruct(k).d) dataStruct(k).d = dataStruct(k).d'; end
    if isrow(dataStruct(k).t) dataStruct(k).t = dataStruct(k).t'; end
    
    % Trim the data. If we include the retract portion of the dataset, then
    % we will go until the applied force detected is negative. If it is not
    % included, we end at the point of maximum z-sensor (maximum
    % deformation condition).
    if ~includeRetract
        dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
        dataStruct(k).d_approach = dataStruct(k).d(1:z_max_ind);
        dataStruct(k).t_approach = dataStruct(k).t(1:z_max_ind);
    else
        F_temp = dataStruct(k).d(z_max_ind:end) .* k_cantilever(k);
        if ~isempty(F_temp) && length(F_temp) > 1
            non_contact_ind = find(F_temp < 0,1);
            if isempty(non_contact_ind) non_contact_ind = length(F_temp)-1; end
        else
            non_contact_ind = 0;
        end
        z_max_ind = non_contact_ind+z_max_ind;
        dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
        dataStruct(k).d_approach = dataStruct(k).d(1:z_max_ind);
        dataStruct(k).t_approach = dataStruct(k).t(1:z_max_ind);
    end
    
    % Give the user feedback, if they weren't aware the file was short.
    if numel(dataStruct(k).z_approach) <= 100
        fprintf('\nWarning: There is a file with very few (<100) z-sensor\ndatapoints in the approach phase:\n%s\n\n',Files(k).name);
    end    
    
    % Calculate Deflection Offset
    [~, d_min_ind] = min(dataStruct(k).d_approach);
    indScale = 0.9;
    d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
    dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
    dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;
    
    % Filter and Shift z_sensor data
    switch filterType
        case 'butter'
            % Create the butterworth
            [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(k).dt)),'low'); % This makes a lowpass filter

            d_smooth = (filter(b,a,dataStruct(k).d_corrected));             % Next, apply the filter
            delay = 0;
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected;
        
        case {'FIR','IIR'}
            Fs = 1/(dataStruct(k).dt); 
            Fstop = ( (round((v_approach(k)),2,'significant')...
                /round(max(v_approach),2,'significant'))...
                /dataStruct(k).dt );
            if Fstop >= 1/(2*dataStruct(k).dt)
                Fstop = 1/(2.05*dataStruct(k).dt);
            elseif Fstop < 1/(10*dataStruct(k).dt)
                Fstop = 1/(10*dataStruct(k).dt);
            end
            Fpass = Fstop*0.01;
            Rp = 0.01;
            Astop = 80;
            LPF = dsp.LowpassFilter('SampleRate',Fs, ...
                                     'FilterType',filterType, ...
                                     'PassbandFrequency',Fpass, ...
                                     'StopbandFrequency',Fstop, ...
                                     'PassbandRipple',Rp, ...
                                     'StopbandAttenuation',Astop);
            delay = floor(mean(grpdelay(LPF)));
        
            d_smooth = LPF(dataStruct(k).d_corrected);

            % Correct filter delay
            sf = d_smooth;
            sf(1:delay) = [];
            
            [~, dSmoothMin] = min(sf);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end

            d_smooth = sf;

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected((delay+1):end);

        case 'none'
            d_smooth = dataStruct(k).d_corrected;
            delay = 0;
            
            [~, dSmoothMin] = min(d_smooth);
        
            if isempty(dSmoothMin) || dSmoothMin <= 0
                dSmoothMin = 1;
            end

            dataStruct(k).z_corrected = dataStruct(k).z_approach - ...
                dataStruct(k).z_approach(dSmoothMin) + ...
                dataStruct(k).d_corrected(dSmoothMin);
            dataStruct(k).dSmoothMin = dSmoothMin;

            z_smooth = dataStruct(k).z_corrected;
        
    end
    
    % Correct the z-sensor data
    dataStruct(k).z_full_corrected = dataStruct(k).z - ...
        dataStruct(k).z_approach(dSmoothMin) + ...
        dataStruct(k).d_corrected(dSmoothMin);

    % Store our detected offsets
    dataStruct(k).d0 = dataStruct(k).d_corrected(dSmoothMin);
    dataStruct(k).z0 = dataStruct(k).z_corrected(dSmoothMin);
    
    % Calculate Force and Indentation
    dataStruct(k).F = dataStruct(k).d_corrected .* k_cantilever(k);
    dataStruct(k).F_smooth = d_smooth .* k_cantilever(k);
    dataStruct(k).d_smooth = d_smooth;
    dataStruct(k).z_smooth = z_smooth;
    dataStruct(k).h = (dataStruct(k).z(1:z_max_ind) - dataStruct(k).z0)...
        - (dataStruct(k).d(1:z_max_ind) - dataStruct(k).d0);
    
    % Get Repulsive Force Application Portion of the Data
    n_offset = length(dataStruct(k).d_corrected(dSmoothMin:z_max_ind));
    n_offset_smooth = length(dataStruct(k).d_corrected(dSmoothMin:(z_max_ind-delay)));
    dt = dataStruct(k).dt;
    
    % Create a "clean" repulsive time array
    t_rep = linspace(0,(n_offset-1)*dt,n_offset)';
    t_rep_smooth = linspace(0,(n_offset_smooth-1)*dt,n_offset_smooth)';

    % Store the repulsive z-sensor and deflection for normal and smooth
    % cases.
    z_rep = dataStruct(k).z_corrected(dSmoothMin:end);
    d_rep = dataStruct(k).d_corrected(dSmoothMin:end);
    z_rep_smooth = dataStruct(k).z_smooth(dSmoothMin:end);
    d_rep_smooth = dataStruct(k).d_smooth(dSmoothMin:end);
    
    % Store the tip position
    tip_rep = d_rep;
    tip_rep_smooth = d_rep_smooth;
    
    % Search for the repulsive portion of the data starting from either the
    % first index and moving to the end ('forward') or starting at the end
    % and walking backward ('reverse').
    if strcmp(findRep,'forward')
        tip_rep_pos = find(tip_rep>0,1);                                   % Find first position above 0
        if isempty(tip_rep_pos) || tip_rep_pos == 0
            tip_rep_pos = 1;
        end
        
        tip_rep_pos_smooth = find(tip_rep_smooth>0,1);                     % Find first position above 0
        if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth == 0
            tip_rep_pos_smooth = 1;
        end
    elseif strcmp(findRep,'reverse')
        tip_rep_pos = (length(tip_rep) - find(flip(tip_rep)<0,1));       % Find last position above 0
        if isempty(tip_rep_pos) || tip_rep_pos == 0
            tip_rep_pos = 1;
        end
        
        tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flip(tip_rep_smooth)<0,1));   % Find last position above 0
        if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth == 0
            tip_rep_pos_smooth = 1;
        end
    end
    
    % Store the repulsive force application portion of the dataset.
    dataStruct(k).t_r = t_rep(tip_rep_pos:end) - t_rep(tip_rep_pos);
    dataStruct(k).z_r = z_rep(tip_rep_pos:end) - z_rep(tip_rep_pos);
    dataStruct(k).d_r = d_rep(tip_rep_pos:end) - d_rep(tip_rep_pos);
    dataStruct(k).t_r_smooth = t_rep_smooth(tip_rep_pos_smooth:end) - t_rep_smooth(tip_rep_pos_smooth);
    dataStruct(k).z_r_smooth = z_rep_smooth(tip_rep_pos_smooth:end) - z_rep_smooth(tip_rep_pos_smooth);
    dataStruct(k).d_r_smooth = d_rep_smooth(tip_rep_pos_smooth:end) - d_rep_smooth(tip_rep_pos_smooth);
    
    dataStruct(k).tip_rep_pos = tip_rep_pos;
    dataStruct(k).tip_rep_pos_smooth = tip_rep_pos_smooth;
    
    % Store for benchmarking later. This will help us decide which region
    % of the curves we can actually "average" together, since we want to
    % use data where all files that are averaged have some contribution
    % (i.e. where they overlap).
    tip_rep_pos_all(k) = tip_rep_pos;
    dSmoothMinAll(k) = dSmoothMin;
    
    % Calculate Force and Indentation during Repulsive Portion
    dataStruct(k).F_r = dataStruct(k).d_r .* k_cantilever(k); % Calculate Force
    dataStruct(k).h_r = (dataStruct(k).z_r - dataStruct(k).d_r); % Calculate Indentation
    
    dataStruct(k).F_r_smooth = dataStruct(k).d_r_smooth .* k_cantilever(k); % Calculate Smooth Force
    dataStruct(k).h_r_smooth = dataStruct(k).z_r_smooth - dataStruct(k).d_r_smooth; % Calculate Smooth Indentation
    
    % Store our critical indices
    dataStruct(k).z_max_ind = z_max_ind;
    dataStruct(k).z_max_ind_smooth = z_max_ind-delay;
    dataStruct(k).dSmoothMin = dSmoothMin;

    % Change to Logarithmic Scaling to use Enrique et al.'s Fit Method
    tr = dataStruct(k).dt;
    st = dataStruct(k).t(z_max_ind);
    
    F_log = [];
    t_log = [];
    
    F_log = log_scale(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
    t_log = log_scale(dataStruct(k).t_r,dataStruct(k).t_r,tr,st);
    
    % Save the Log-Form of the tip force (F) and time array (t_r)
    dataStruct(k).F_r_log = F_log;
    dataStruct(k).t_r_log = t_log;

end

clearvars k point_number_temp run_number_temp

% Determine how many load levels we are considering
v_approach = round(v_approach,2,'significant');
v_unique = (unique(v_approach));

if length(Files) > 1
    N = length(Files);
    
    for k=1:length(Files)
        [minVal(k),~] = min(dataStruct(k).t_approach);
        [maxVal(k),~] = max(dataStruct(k).t_approach);
        dtVal(k) = dataStruct(k).dt;
    end
    
    if length(v_unique) == 1
        
        % Find files corresponding to current velocity
        [startNum,~] = max(minVal);
        [minRepInd,minRepFile] = min(tip_rep_pos_all + dSmoothMinAll);
        [endNum,~] = min(maxVal);

        xi = startNum:median(dtVal):endNum;  % Create Vector Of Common time values
        di = [];
        zi = [];
        ti = [];
        dataLengths = zeros(1,length(Files));
        currentFile = zeros(1,length(Files));
        forceLimits = zeros(1,length(Files));
        timeLimits = zeros(1,length(Files));
                
        % When only one velocity is present
        for k=1:length(Files)
            shiftVal = (dataStruct(k).t_approach(tip_rep_pos_all(k)+dataStruct(k).dSmoothMin) - ...
                            dataStruct(minRepFile).t_approach(minRepInd));

            [tempz, ~] = unique(dataStruct(k).t_approach - shiftVal);

            di(k,:) = interp1(tempz, dataStruct(k).d_corrected,...
                xi(:), 'linear', NaN); % Interploate deflection to new �x� Values

            zi(k,:) = interp1(tempz, dataStruct(k).z_corrected,...
                xi(:), 'linear', NaN); % Interploate z-sensor to new �x� Values

            % Create an associated time array for averaging
            ti(k,:) =  xi;
            
            % Hold on to the data lengths so we can keep track
            % of which files are worst and should be ignored.
            if isempty(find(isnan(zi(k,:)),1,'first'))
                dataLengths(k) = size(zi(k,:),2);
                forceLimits(k) = max(di(k,:));
                timeLimits(k) = max(tempz);
            else
                dataLengths(k) = find(isnan(zi(k,:)),1,'first');
                forceLimits(k) = max(di(k,:));
                timeLimits(k) = max(tempz);
            end
            currentFile(k) = k;
            
        end
        
        % Interpolate according to the variance in the time array that
        % comes from averaging the results.
        t_interp = (1:size(ti,2))*median(dtVal);
        z_interp = interp1(mean(ti,1), mean(zi,1),...
                t_interp, 'linear', NaN)'; % Interploate Or Extrapolate To New Time Values
        d_interp = interp1(mean(ti,1), mean(di,1),...
                t_interp, 'linear', NaN)'; % Interploate Or Extrapolate To New Time Values
        
        nanCheck = isnan(reshape(t_interp,1,length(t_interp))) +...
            isnan(reshape(z_interp,1,length(z_interp))) +...
            isnan(reshape(d_interp,1,length(d_interp)));
        nanCheck(nanCheck ~= 0) = 1;

        t_interp(nanCheck == 1) = [];
        z_interp(nanCheck == 1) = [];
        d_interp(nanCheck == 1) = [];
        
        if size(z_interp,2) < size(z_interp,1)
            dataStruct(length(Files)+1).z_average = z_interp';
        else
            dataStruct(length(Files)+1).z_average = z_interp;
        end
        
        if size(d_interp,2) < size(d_interp,1)
            dataStruct(length(Files)+1).d_average = d_interp';
        else
            dataStruct(length(Files)+1).d_average = d_interp;
        end
        
        if size(t_interp,2) < size(t_interp,1)
            dataStruct(length(Files)+1).t_average = t_interp';
        else
            dataStruct(length(Files)+1).t_average = t_interp;
        end

        % Sanity check --- should be the same as the datasets' dt
        dataStruct(length(Files)+1).dt = mean(diff(dataStruct(length(Files)+1).t_average));

    else
                
        % Loop through unique velocities
        for j=1:length(v_unique)
            
            % Find files corresponding to current velocity
            velInd = zeros(1,N);
            velInd(v_approach == v_unique(j)) = 1;
            
            [startNum,~] = max(minVal(velInd==1));
            [minRepInd,temp] = min(tip_rep_pos_all(velInd==1) + dSmoothMinAll(velInd==1));
            fileInds = find(velInd);
            minRepFile = fileInds(temp);
            [endNum,~] = min(maxVal(velInd==1));
                        
            xi = (startNum:median(dtVal):endNum)';  % Create Vector Of Common time values
            di = [];
            zi = [];
            ti = [];
            dataLengths = zeros(1,length(Files));
            currentFile = zeros(1,length(Files));
            forceLimits = zeros(1,length(Files));
            timeLimits = zeros(1,length(Files));
                        
            if sum(v_approach(:) == v_unique(j)) > 1
                % Loop through all files
                for k=1:length(Files)
                    % Grab the tip radius for this velocity and 

                    % Check if this file is relevant for this specific
                    % averaging operation
                    if v_approach(k) == v_unique(j)
                        shiftVal = (dataStruct(k).t_approach(tip_rep_pos_all(k)+dataStruct(k).dSmoothMin) - ...
                            dataStruct(minRepFile).t_approach(minRepInd));
                        
                        [tempz, ~] = unique(dataStruct(k).t_approach - shiftVal);

                        di = horzcat(di, interp1(tempz, dataStruct(k).d_corrected,...
                            xi(:), 'linear', NaN)); % Interploate deflection to new �x� Values

                        zi = horzcat(zi, interp1(tempz, dataStruct(k).z_corrected,...
                            xi(:), 'linear', NaN)); % Interploate z-sensor to new �x� Values

                        % Create an associated time array for averaging
                        ti =  horzcat(ti, xi);
                        
                        % Hold on to the data lengths so we can keep track
                        % of which files are worst and should be ignored.
                        if isempty(find(isnan(zi(:,end)),1,'first'))
                            dataLengths = horzcat(dataLengths, size(zi(:,end),2));
                            forceLimits = horzcat(forceLimits, max(di(:,end)));
                            timeLimits = horzcat(timeLimits, max(tempz));
                        else
                            dataLengths = horzcat(dataLengths, find(isnan(zi(:,end)),1,'first'));
                            forceLimits = horzcat(forceLimits, max(di(:,end)));
                            timeLimits = horzcat(timeLimits, max(tempz));
                        end
                        currentFile = horzcat(currentFile, k);
                                                                        
                    else
                        continue;
                    end

                end
            else
                % There is only one file for this velocity, so treat it
                % like a single file would be then skip to the next
                % iteration.
                velIndRef = find(velInd);
                
                if size(dataStruct(velIndRef).z_corrected,2) < size(dataStruct(velIndRef).z_corrected,1)
                    dataStruct(length(Files)+j).z_average = dataStruct(velIndRef).z_corrected';
                else
                    dataStruct(length(Files)+j).z_average = dataStruct(velIndRef).z_corrected;
                end
                
                if size(dataStruct(velIndRef).d_corrected,2) < size(dataStruct(velIndRef).d_corrected,1)
                    dataStruct(length(Files)+j).d_average = dataStruct(velIndRef).d_corrected';
                else
                    dataStruct(length(Files)+j).d_average = dataStruct(velIndRef).d_corrected;
                end
                
                if size(dataStruct(velIndRef).t_approach,2) < size(dataStruct(velIndRef).t_approach,1)
                    dataStruct(length(Files)+j).t_average = dataStruct(velIndRef).t_approach';
                else
                    dataStruct(length(Files)+j).t_average = dataStruct(velIndRef).t_approach;
                end
                
                dataStruct(length(Files)+j).dt = dataStruct(velIndRef).dt;
                continue;
                
            end
            
            % Average this load level's curves.
            % Interpolate according to the variance in the time array that
            % comes from averaging the results.
            t_interp = (1:size(ti,1))*median(dtVal)';
            z_interp = interp1(mean(ti,2), mean(zi,2),...
                    t_interp, 'linear', NaN); % Interploate To New Time Values
            d_interp = interp1(mean(ti,2), mean(di,2),...
                    t_interp, 'linear', NaN); % Interploate To New Time Values
            
            nanCheck = isnan(reshape(t_interp,length(t_interp),1)) +...
                isnan(reshape(z_interp,length(z_interp),1)) +...
                isnan(reshape(d_interp,length(d_interp),1));
            nanCheck(nanCheck ~= 0) = 1;
            
            t_interp(nanCheck == 1) = [];
            z_interp(nanCheck == 1) = [];
            d_interp(nanCheck == 1) = [];
                
            if size(z_interp,2) < size(z_interp,1)
                dataStruct(length(Files)+j).z_average = z_interp';
            else
                dataStruct(length(Files)+j).z_average = z_interp;
            end

            if size(d_interp,2) < size(d_interp,1)
                dataStruct(length(Files)+j).d_average = d_interp';
            else
                dataStruct(length(Files)+j).d_average = d_interp;
            end

            if size(t_interp,2) < size(t_interp,1)
                dataStruct(length(Files)+j).t_average = t_interp';
            else
                dataStruct(length(Files)+j).t_average = t_interp;
            end
            
            % Sanity check --- should be the same as the datasets' dt
            dataStruct(length(Files)+j).dt = mean(diff(dataStruct(length(Files)+1).t_average));
            
        end
    end
    
else
    % There is only one file being analyzed, so simply copy over the data
    % to the "average" row.
    dataStruct(length(Files)+1).z_average = dataStruct(1).z_corrected;
    dataStruct(length(Files)+1).d_average = dataStruct(1).d_corrected;
    dataStruct(length(Files)+1).t_average = dataStruct(1).t_approach;
    dataStruct(length(Files)+1).dt = dataStruct(1).dt;
end

% Pre-Processing for Averaged Data of all load levels.
for i = 1:length(v_unique)
        
    % Find the approach portion of the data. This is only really necessary
    % if the average data happens to go past the limit of all the datasets
    % for some reason.
    [~, z_max_ind] = max(dataStruct(length(Files)+i).z_average);
    
    if ~includeRetract
        dataStruct(length(Files)+i).z_approach = dataStruct(length(Files)+i).z_average(1:z_max_ind);
        dataStruct(length(Files)+i).d_approach = dataStruct(length(Files)+i).d_average(1:z_max_ind);
    else
        F_temp = dataStruct(length(Files)+i).d_average(z_max_ind:end) .* mean(k_cantilever);
        if ~isempty(F_temp) && length(F_temp) > 1
            non_contact_ind = find(F_temp < 0,1);
            if isempty(non_contact_ind) non_contact_ind = length(F_temp)-1; end
        else
            non_contact_ind = 0;
        end
        z_max_ind = non_contact_ind+z_max_ind;
        dataStruct(length(Files)+i).z_approach = dataStruct(length(Files)+i).z_average(1:z_max_ind);
        dataStruct(length(Files)+i).d_approach = dataStruct(length(Files)+i).d_average(1:z_max_ind);
    end
    
    % Filter and Shift z_sensor data
    switch filterType
        case 'butter'
            % Create the butterworth
            [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(length(Files)+i).dt)),'low'); % This makes a lowpass filter
            d_approach_smooth = (filter(b,a,dataStruct(length(Files)+i).d_approach)); % Next, apply the filter
    
            [~, dSmoothMin] = min(d_approach_smooth);   
            z_approach_smooth = dataStruct(length(Files)+i).z_approach;
        
        case {'IIR','FIR'}
            % Use an IIR or FIR filter on the data
            Fs = 1/(dataStruct(length(Files)+i).dt); 
            Fstop = 100*( (round((v_unique(i)),2,'significant')...
                /round(max(v_unique),2,'significant'))...
                /dataStruct(length(Files)+i).dt );
            if Fstop >= 1/(3*dataStruct(length(Files)+i).dt)
                Fstop = 1/(3*dataStruct(length(Files)+i).dt);
            elseif Fstop < 1/(10*dataStruct(length(Files)+i).dt)
                Fstop = 1/(10*dataStruct(length(Files)+i).dt);
            end
            Fpass = Fstop*0.01;
            Rp = 0.01;
            Astop = 80;
            LPF = dsp.LowpassFilter('SampleRate',Fs, ...
                                     'FilterType',filterType, ...
                                     'PassbandFrequency',Fpass, ...
                                     'StopbandFrequency',Fstop, ...
                                     'PassbandRipple',Rp, ...
                                     'StopbandAttenuation',Astop);
            delay = floor(mean(grpdelay(LPF)));

            d_approach_smooth = LPF(dataStruct(length(Files)+i).d_approach);         % Smooth with IIR filter

            % Correct filter delay
            sf = d_approach_smooth;
            sf(1:delay) = [];
            
            [~, dSmoothMin] = min(sf);
            d_approach_smooth = sf;
            z_approach_smooth = dataStruct(length(Files)+i).z_approach((delay+1):end);

        case 'none'
            % Apply no filtering
            d_approach_smooth = dataStruct(length(Files)+i).d_approach;
            delay = 0;
            
            [~, dSmoothMin] = min(d_approach_smooth);   
            z_approach_smooth = dataStruct(length(Files)+i).z_approach;
        
    end
       
    t_rep = dataStruct(length(Files)+i).t_average(dSmoothMin:z_max_ind);
    z_rep = dataStruct(length(Files)+i).z_average(dSmoothMin:z_max_ind);
    d_rep = dataStruct(length(Files)+i).d_average(dSmoothMin:z_max_ind);
    t_rep_smooth = dataStruct(length(Files)+i).t_average(dSmoothMin:(z_max_ind-delay));
    z_rep_smooth = z_approach_smooth(dSmoothMin:(z_max_ind-delay));
    d_rep_smooth = d_approach_smooth(dSmoothMin:(z_max_ind-delay));
    
    tip_rep = d_rep;
    tip_rep_smooth = d_rep_smooth;
    if strcmp(findRep,'forward')
        tip_rep_pos = find(tip_rep>0,1);                                   % Find first position above 0
        tip_rep_pos_smooth = find(tip_rep_smooth>0,1);                     % Find first position above 0
    elseif strcmp(findRep,'reverse')
        tip_rep_pos = (length(tip_rep) - find(flip(tip_rep)<0,1));       % Find last position above 0
        if isempty(tip_rep_pos) || tip_rep_pos == 0
            tip_rep_pos = 1;
        end
        
        tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flip(tip_rep_smooth)<0,1));       % Find last position above 0
        if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth == 0
            tip_rep_pos_smooth = 1;
        end
    end
    
    % Find the repulsive portion (force application) region
    dataStruct(length(Files)+i).t_r = t_rep(tip_rep_pos:end) - t_rep(tip_rep_pos);
    dataStruct(length(Files)+i).z_r = z_rep(tip_rep_pos:end) - z_rep(tip_rep_pos);
    dataStruct(length(Files)+i).d_r = d_rep(tip_rep_pos:end) - d_rep(tip_rep_pos);
    t_r_smooth = t_rep_smooth(tip_rep_pos_smooth:end) - t_rep_smooth(tip_rep_pos_smooth);
    z_r_smooth = z_rep_smooth(tip_rep_pos_smooth:end) - z_rep_smooth(tip_rep_pos_smooth);
    d_r_smooth = d_rep_smooth(tip_rep_pos_smooth:end) - d_rep_smooth(tip_rep_pos_smooth);
    
    % Calculate Force and Indentation during Repulsive Portion
    dataStruct(length(Files)+i).F_r = dataStruct(length(Files)+i).d_r .* mean(k_cantilever);
    dataStruct(length(Files)+i).h_r = dataStruct(length(Files)+i).z_r - dataStruct(length(Files)+i).d_r; % Calculate Indentation
    
    dataStruct(length(Files)+i).z_max_ind = z_max_ind;
    dataStruct(length(Files)+i).z_max_ind_smooth = z_max_ind - delay;
    dataStruct(length(Files)+i).dSmoothMin = dSmoothMin;
    
    % Create smooth Force and Indentation
    k = 1;
    k_avg = 0;
    while k_avg == 0
        % Check if this file is relevant for finding cantilever stiffness
        if v_approach(k) == v_unique(i)
            k_avg = k_cantilever(i);
            continue;
        else
            k = k+1;
            continue;
        end
    end
    
    dataStruct(length(Files)+i).t_r_smooth = t_r_smooth;
    dataStruct(length(Files)+i).z_r_smooth = z_r_smooth;
    dataStruct(length(Files)+i).d_r_smooth = d_r_smooth;
    dataStruct(length(Files)+i).F_r_smooth = d_r_smooth .* k_avg; % Calculate Smooth Force
    dataStruct(length(Files)+i).h_r_smooth = z_r_smooth - d_r_smooth; % Calculate Smooth Indentation
    
    if removeNegatives
        % Original
        toRemove = (dataStruct(length(Files)+i).h_r <= 0 | dataStruct(length(Files)+i).F_r <= 0);
        dataStruct(length(Files)+i).t_r(toRemove) = [];
        dataStruct(length(Files)+i).z_r(toRemove) = [];
        dataStruct(length(Files)+i).d_r(toRemove) = [];
        dataStruct(length(Files)+i).F_r(toRemove) = [];
        dataStruct(length(Files)+i).h_r(toRemove) = [];
    
        % Smooth
        toRemove = (dataStruct(length(Files)+i).h_r_smooth <= 0 | dataStruct(length(Files)+i).F_r_smooth <= 0);
        dataStruct(length(Files)+i).t_r_smooth(toRemove) = [];
        dataStruct(length(Files)+i).z_r_smooth(toRemove) = [];
        dataStruct(length(Files)+i).d_r_smooth(toRemove) = [];
        dataStruct(length(Files)+i).F_r_smooth(toRemove) = [];
        dataStruct(length(Files)+i).h_r_smooth(toRemove) = [];
    end
    
    % Change to Logarithmic Scaling to use Enrique et al.'s Fit Method
    tr = dataStruct(length(Files)+i).dt;
    st = dataStruct(length(Files)+i).t_r(end);

    F_log = [];
    t_log = [];

    F_log = log_scale(dataStruct(length(Files)+i).F_r,dataStruct(length(Files)+i).t_r,tr,st);
    t_log = log_scale(dataStruct(length(Files)+i).t_r,dataStruct(length(Files)+i).t_r,tr,st);
    
    % Save the Log-Form of the tip force (F) and time array (t_r)
    dataStruct(length(Files)+i).F_r_log = F_log;
    dataStruct(length(Files)+i).t_r_log = t_log;

end

end

