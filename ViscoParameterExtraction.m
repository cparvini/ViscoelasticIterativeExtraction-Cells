format long
clear all
close all
clc

% ======================== %
% Viscoelastic Parameter
% Extraction Script
%
% Created by cparvini
% ======================== %
% This file will perform the 
% required data analysis on 
% AFM SFS curves to extract 
% meaningful information about 
% how stiff the samples are on 
% specific timescales.
% ======================== %

addpath(genpath([pwd '/includes']))
addpath data

%% Settings for Data Analysis
% Get Experiment Settings 
r_tip = 12.5e-6;                        % DOUBLE, Tip Radius in METERS
nu_sample = 0.5;                        % DOUBLE, Poisson's ratio of sample. For incompressible, set to 0.5

% Filter and Data Correction settings
removeNegatives = 1;                    % (1) Remove negative values from the input data (due to noise); (0) Do not remove negative values
filterType = 'FIR';                     % ('none') Save the unfiltered data into the "smooth data" fields (not helpful); ('FIR') Use an FIR filter on the data; ('butter') Use an Nth order butterworth on the data
findRep = 'forward';                    % ('forward') Find the repulsive portion of the dataset starting from t=0; ('reverse') Find the repulsive portion of the dataset starting from t=t(end)
N = 2;                                  % INTEGER, Filter order for the Butterworth (if using, ignored if not)
cutoff_Hz = 5e3;                        % DOUBLE, Butterworth cutoff frequency (if using, ignored if not)

% Viscoelastic Model Settings
n_terms_tot = 3;                        % Maximum number of terms to allow in the viscoelastic series
elasticSetting = 'y';                   % ('y') include an elastic element in the generalized voigt model; ('n') exclude the element
elasticSetting_maxwell = 'y';           % ('y') include an elastic element in the generalized maxwell model; ('n') exclude the element
fluidSetting = 'n';                     % ('y') include a fluidity element in the generalized voigt model; ('n') exclude the element
fluidSetting_maxwell = 'n';             % ('y') include a fluidity element in the generalized maxwell model; ('n') exclude the element

% Manual Simulation Setting
useBruker = 0;                          % (1) Use the Bruker toolbox to read your .spm files (MUST HAVE TOOLBOX INSTALLED FIRST); (0) Don't use the Bruker toolbox
n_samples = 100;                        % INTEGER, The number of random-starts to attempt per model configuration
n_maxIterations = 1e4;                  % INTEGER, The number of internal algorithm iterations allowed per random-start
multiSim = 'y';                         % ('y') Run the simulation seperately for terms 1-N where N==n_terms_tot; ('n') Only run for N==n_terms_tot
skipNums = 0;                           % INTEGER ARRAY; each member of the array is a model configuration that will be SKIPPED. Only useful if multiSim=='y'        
avAnalysis = 'a';                       % ('a') Use the AVERAGE of all force curves with the SAME APPROACH VELOCITY for fitting; ('i') Process/fit files INDIVIDUALLY
smoothData = 'y';                       % ('y') Use the FILTERED data for analysis; ('n') Use the UNFILTERED data for analysis
scaling = 1e-60;                        % DOUBLE, Convergence tolerance for the optimization algorithm
forwardFitTimescale = 1;                % (1) Introduce characteristic times from SHORTER to LONGER times; (2) Introduce characteristic times from LONGER to SHORTER times (based on length of experiment)
minTimescale = 1e-4;                    % DOUBLE, Minimum timescale to use for the first element
timeOperation = 0;                      % (1) Save the time required for different steps, to evaluate performance later; (2) Do not record time to monitor performance
fitMethod = 'iterative';                % ('iterative') Use iterative term introduction (preferred); ('open') Use the open-search methodology (LEGACY)
fitLog = 1;                             % (1) Fit using Log-Sampled data; (0) Do NOT use Log-Sampled Data for Fitting
tossThresh = 0;                         % DOUBLE, Throw away (set to NaN) parameter sets that are within (tossThresh) orders of magnitude of the upper/lower bounds
clipData = 1;                           % (1) During iterative term introduction, use only data that can be represented by the characteristic times currently within the model configuration; (0) Use the full-fidelity time data for ALL configurations
includeRetract = 0;                     % (1) Include the retract data in the fitting process (NOT RECOMMENDED, it is hard to determine the point where contact is lost); (0) Exclude retract data from the fitting process
if (includeRetract) clipData = 0; end   % You can't use the retract if you clip data.
fitElasticFirst = 1;                    % (1) During iterative term introduction, fit the elastic term separately (if included above); (0) Do not fit the elastic term separately---instead, fit it at the same time as the first element

% Save the settings for later
settingsStruct = struct;
settingsStruct.removeNegatives = removeNegatives;
settingsStruct.filterType = filterType;
settingsStruct.findRep = findRep;
settingsStruct.N = N;
settingsStruct.cutoff_Hz = cutoff_Hz;
settingsStruct.r_tip = r_tip;
settingsStruct.useBruker = useBruker;
settingsStruct.n_terms_tot = n_terms_tot;
settingsStruct.n_samples = n_samples;
settingsStruct.n_maxIterations = n_maxIterations;
settingsStruct.multiSim = multiSim;
settingsStruct.skipNums = skipNums;
settingsStruct.elasticSetting = elasticSetting;
settingsStruct.elasticSetting_maxwell = elasticSetting_maxwell;
settingsStruct.fluidSetting = fluidSetting;
settingsStruct.fluidSetting_maxwell = fluidSetting_maxwell;
settingsStruct.avAnalysis = avAnalysis;
settingsStruct.smoothData = smoothData;
settingsStruct.scaling = scaling;
settingsStruct.nu_sample = nu_sample;
settingsStruct.forwardFitTimescale = forwardFitTimescale;
settingsStruct.minTimescale = minTimescale;
settingsStruct.timeOperation = timeOperation;
settingsStruct.fitMethod = fitMethod;
settingsStruct.clipData = clipData;
settingsStruct.fitElasticFirst = fitElasticFirst;

%% Load AFM Files and Parse
% Ask the user to select their experiment data directory
pathname = uigetdir(pwd,...
    'Select the Folder Containing your AFM Files');
savePrepend = input('Please Enter the Label for the Output Files: ','s');

% Save our settings file
save([pathname sprintf('/%s-settingsStruct-average.mat',savePrepend)],'settingsStruct')
fprintf('\nSaved the simulation settings to an m-file.\n');
clearvars settingsStruct

% Perform the data loading steps
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
        error('ERROR: Attempt to pass multiple experiments into the script has FAILED. Please ensure only one experiment label is used in this directory.');
    end
else
    error('ERROR: No experiment files found! Please ensure you have selected the correct directory.');
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
            run_number = str2num(strrep(FileInfo{4},'-','.'))+1;         % Run Number
            
    end
end

% uniqueVelocities = unique(v_approach);
dataStruct = struct;

% Remove the files we don't care about
FilesRemove=(~ismember({Files.name},{FilesCheck.name}));
Files(FilesRemove) = [];

for k=1:length(Files)
    FileNames = Files(k).name;
    FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
    
    strSwitch = FileInfo{1};
    while isstrprop(strSwitch(end),'digit')
        strSwitch(end) = [];
    end

    switch lower(strSwitch)
        case lower('FD')
            RawData = IBWread([Files(k).folder '/' Files(k).name]);
            [headerValue,null] = strsplit(RawData.WaveNotes,'\r',...
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
            
            if useBruker == 1
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
                
                [time, deflectionError, xLabel, yLabel] = NSMU.CreateForceTimePlot(1, NSMU.VOLTS);

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
                
                nu_sample = NSMU.GetPoissonRatio();
                r_tip = NSMU.GetTipRadius();
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
                
                approachcurvesraw=flipud(usefulData(a_ind(1):b_ind(1)));
                retractcurvesraw=flipud(usefulData(a_ind(2):b_ind(2)));
                zapproachraw=flipud(usefulData(a_ind(3):b_ind(3)));
                zretractraw=flipud(usefulData(a_ind(4):b_ind(4)));

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
                dataError_pos = ( length(zApproach) - find(flipud(gradient(zApproach(1:end-buffer)))<0,1) );    % Find last position above 0
                dataError_pos2 = ( length(approachDefl) - find(abs(flipud(gradient(approachDefl)))>10*max(abs(gradient(approachDefl(end-buffer:end)))),1) );    % Find where deflection blows up
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
                    dataStruct(k).t = time(1:(end-dataError_pos+buffer));
                    dataStruct(k).dt = dt;
                    dataStruct(k).n_sam = length(dataStruct(k).t);
                end
                
                v_approach(k) = mean(abs(diff(dataStruct(k).z))./(dataStruct(k).dt*10));    % Approach Velocity, m/s
                fclose(fid);
                
            end
            
    end
    
    % Pre-Processing
    % Find the approach portion of the data    
    [z_max, z_max_ind] = max(dataStruct(k).z);
    
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

    if size(dataStruct(k).z_approach,1) <= 100
        fprintf('\nThere is a bad file, with very few z-sensor datapoints:\n%s\n\n',Files(k).name);
    end    
    
    % Calculate Deflection Offset
    [null, d_min_ind] = min(dataStruct(k).d_approach);
    indScale = 0.9;
    d_0_mean = mean(dataStruct(k).d_approach(1:round(d_min_ind*indScale)));
    dataStruct(k).d_corrected = dataStruct(k).d_approach - d_0_mean;
    dataStruct(k).d_full_corrected = dataStruct(k).d - d_0_mean;
    
    % Filter and Shift z_sensor data
    if strcmp(filterType,'butter')
        % Create the butterworth
        [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(k).dt)),'low'); % This makes a lowpass filter

        d_smooth = (filter(b,a,dataStruct(k).d_corrected));             % Next, apply the filter
        delay = 0;
        
    elseif strcmp(filterType,'FIR')
        
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
        
    elseif strcmp(filterType,'movingAverageFFT')
        
        % Filtering using a moving average in the frequency domain
        f = linspace(0,1/(dataStruct(k).dt),length(dataStruct(k).d_corrected));
        data = (fftshift(fft(dataStruct(k).d_corrected)));
        M = movmean(data,7,'SamplePoints',f);
        
        figure
        ftBerk(dataStruct(k).t_approach,dataStruct(k).d_corrected,1);
        hold on
        ftBerk(dataStruct(k).t_approach,abs(ifft(ifftshift(M))),1);
        hold off
        
        d_smooth = abs(ifft(ifftshift(M)));
        delay = 0;
        
        % Force findrep to move in reverse...
        % Necessary or it won't find the right minimum
        findRep = 'reverse';
        
    elseif strcmp(filterType,'none')
        
        d_smooth = dataStruct(k).d_corrected;
        delay = 0;
        
    end
        
    if exist('sf','var')
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
    else
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
    
    dataStruct(k).z_full_corrected = dataStruct(k).z - ...
        dataStruct(k).z_approach(dSmoothMin) + ...
        dataStruct(k).d_corrected(dSmoothMin);

    dataStruct(k).d0 = dataStruct(k).d_corrected(dSmoothMin);
    dataStruct(k).z0 = dataStruct(k).z_corrected(dSmoothMin);
    
    % Calculate Force and Indentation
    dataStruct(k).F = dataStruct(k).d_corrected .* k_cantilever(k);
    dataStruct(k).F_smooth = d_smooth .* k_cantilever(k);
    dataStruct(k).d_smooth = d_smooth;
    dataStruct(k).z_smooth = z_smooth;
    dataStruct(k).h = (dataStruct(k).z(1:z_max_ind) - dataStruct(k).z0)...
        - (dataStruct(k).d(1:z_max_ind) - dataStruct(k).d0); % Calculate Indentation
    
    % Get Repulsive Portion of the Data
    n_offset = length(dataStruct(k).d_corrected(dSmoothMin:z_max_ind));
    n_offset_smooth = length(dataStruct(k).d_corrected(dSmoothMin:(z_max_ind-delay)));
    dt = dataStruct(k).dt;
    
    t_rep = linspace(0,(n_offset-1)*dt,n_offset)';
    t_rep_smooth = linspace(0,(n_offset_smooth-1)*dt,n_offset_smooth)';

    z_rep = dataStruct(k).z_corrected(dSmoothMin:end);
    d_rep = dataStruct(k).d_corrected(dSmoothMin:end);
    z_rep_smooth = dataStruct(k).z_smooth(dSmoothMin:end);
    d_rep_smooth = dataStruct(k).d_smooth(dSmoothMin:end);
    
    tip_rep = d_rep;
    tip_rep_smooth = d_rep_smooth;
    
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
        tip_rep_pos = (length(tip_rep) - find(flipud(tip_rep)<0,1));       % Find last position above 0
        if isempty(tip_rep_pos) || tip_rep_pos == 0
            tip_rep_pos = 1;
        end
        
        tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flipud(tip_rep_smooth)<0,1));   % Find last position above 0
        if isempty(tip_rep_pos_smooth) || tip_rep_pos_smooth == 0
            tip_rep_pos_smooth = 1;
        end
    end
    
    dataStruct(k).t_r = t_rep(tip_rep_pos:end) - t_rep(tip_rep_pos);
    dataStruct(k).z_r = z_rep(tip_rep_pos:end) - z_rep(tip_rep_pos);
    dataStruct(k).d_r = d_rep(tip_rep_pos:end) - d_rep(tip_rep_pos);
    dataStruct(k).t_r_smooth = t_rep_smooth(tip_rep_pos_smooth:end) - t_rep_smooth(tip_rep_pos_smooth);
    dataStruct(k).z_r_smooth = z_rep_smooth(tip_rep_pos_smooth:end) - z_rep_smooth(tip_rep_pos_smooth);
    dataStruct(k).d_r_smooth = d_rep_smooth(tip_rep_pos_smooth:end) - d_rep_smooth(tip_rep_pos_smooth);
    
    dataStruct(k).tip_rep_pos = tip_rep_pos;
    dataStruct(k).tip_rep_pos_smooth = tip_rep_pos_smooth;
    tip_rep_pos_all(k) = tip_rep_pos;
    tip_rep_pos_all_smooth(k) = tip_rep_pos_smooth;
    dSmoothMinAll(k) = dSmoothMin;
    
    % Calculate Force and Indentation during Repulsive Portion
    dataStruct(k).F_r = dataStruct(k).d_r .* k_cantilever(k); % Calculate Force
    dataStruct(k).h_r = (dataStruct(k).z_r - dataStruct(k).d_r); % Calculate Indentation
    
    dataStruct(k).F_r_smooth = dataStruct(k).d_r_smooth .* k_cantilever(k); % Calculate Smooth Force
    dataStruct(k).h_r_smooth = dataStruct(k).z_r_smooth - dataStruct(k).d_r_smooth; % Calculate Smooth Indentation
    
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

% Determine how many load levels we are considering
v_raw = v_approach;
v_approach = round(v_approach,2,'significant');
v_unique = (unique(v_approach));
[~,indv] = sort(v_approach);

if length(Files) > 1
    FileNames = Files(1).name;
    FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

    N = length(Files);
    
    for k=1:length(Files)
        [minVal(k), indMin(k)] = min(dataStruct(k).t_approach);
        [maxVal(k), indMax(k)] = max(dataStruct(k).t_approach);
        numPoints(k) = length(dataStruct(k).t_approach(indMin(k):indMax(k)));
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
        dataLengths = [];
        currentFile = [];
        forceLimits = [];
        timeLimits = [];
        
        % When only one velocity is present
        for k=1:length(Files)
            shiftVal = (dataStruct(k).t_approach(tip_rep_pos_all(k)+dataStruct(k).dSmoothMin) - ...
                            dataStruct(minRepFile).t_approach(minRepInd));

            [tempz, indz] = unique(dataStruct(k).t_approach - shiftVal);

            di(:,k) = interp1(tempz, dataStruct(k).d_corrected,...
                xi(:), 'linear', NaN); % Interploate deflection to new �x� Values

            zi(:,k) = interp1(tempz, dataStruct(k).z_corrected,...
                xi(:), 'linear', NaN); % Interploate z-sensor to new �x� Values

            % Create an associated time array for averaging
            ti(:,k) =  xi';
            
            % Hold on to the data lengths so we can keep track
            % of which files are worst and should be ignored.
            if isempty(find(isnan(zi(:,k)),1,'first'))
                dataLengths(k) = size(zi(:,k),1);
                forceLimits(k) = max(di(:,k));
                timeLimits(k) = max(tempz);
            else
                dataLengths(k) = find(isnan(zi(:,k)),1,'first');
                forceLimits(k) = max(di(:,k));
                timeLimits(k) = max(tempz);
            end
            currentFile(k) = k;
            
        end
        
        % Keep track of the shortest file
        [minLength,minInd] = min(dataLengths);
        [~,sortInd] = sort(dataLengths);
        [~,sortInd2] = sort(forceLimits);
        [~,sortInd3] = sort(timeLimits);
        dataStruct(length(Files)+1).shortestLength = minLength;
        dataStruct(length(Files)+1).shortestFile = currentFile(minInd);
        dataStruct(length(Files)+1).sortedFiles = currentFile(sortInd);
        dataStruct(length(Files)+1).sortedLengths = dataLengths(sortInd);
        dataStruct(length(Files)+1).sortedForceLimits = forceLimits(sortInd2);
        dataStruct(length(Files)+1).sortedForceFiles = currentFile(sortInd2);
        dataStruct(length(Files)+1).sortedTimeLimits = timeLimits(sortInd3);
        dataStruct(length(Files)+1).sortedTimeFiles = currentFile(sortInd3);

        % Interpolate according to the variance in the time array that
        % comes from averaging the results.
        t_interp = (1:size(ti,1))*median(dtVal);
        z_interp = interp1(mean(ti,2), mean(zi,2),...
                t_interp, 'linear', NaN)'; % Interploate Or Extrapolate To New Time Values
        d_interp = interp1(mean(ti,2), mean(di,2),...
                t_interp, 'linear', NaN)'; % Interploate Or Extrapolate To New Time Values
        
        nanCheck = isnan(reshape(t_interp,length(t_interp),1)) +...
            isnan(reshape(z_interp,length(z_interp),1)) +...
            isnan(reshape(d_interp,length(d_interp),1));
        nanCheck(nanCheck ~= 0) = 1;

        t_interp(nanCheck == 1) = [];
        z_interp(nanCheck == 1) = [];
        d_interp(nanCheck == 1) = [];
            
        if size(z_interp,2) > size(z_interp,1)
            dataStruct(length(Files)+1).z_average = z_interp';
        else
            dataStruct(length(Files)+1).z_average = z_interp;
        end
        
        if size(d_interp,2) > size(d_interp,1)
            dataStruct(length(Files)+1).d_average = d_interp';
        else
            dataStruct(length(Files)+1).d_average = d_interp;
        end
        
        if size(t_interp,2) > size(t_interp,1)
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
                        
            xi = startNum:median(dtVal):endNum;  % Create Vector Of Common time values
            di = [];
            zi = [];
            ti = [];
            dataLengths = [];
            currentFile = [];
            forceLimits = [];
            timeLimits = [];
            
            if sum(v_approach(:) == v_unique(j)) > 1
                % Loop through all files
                for k=1:length(Files)

                    % Check if this file is relevant for this specific
                    % averaging operation
                    if v_approach(k) == v_unique(j)
                        shiftVal = (dataStruct(k).t_approach(tip_rep_pos_all(k)+dataStruct(k).dSmoothMin) - ...
                            dataStruct(minRepFile).t_approach(minRepInd));
                        
                        [tempz, indz] = unique(dataStruct(k).t_approach - shiftVal);

                        di = horzcat(di, interp1(tempz, dataStruct(k).d_corrected,...
                            xi(:), 'linear', NaN)); % Interploate deflection to new �x� Values

                        zi = horzcat(zi, interp1(tempz, dataStruct(k).z_corrected,...
                            xi(:), 'linear', NaN)); % Interploate z-sensor to new �x� Values

                        % Create an associated time array for averaging
                        ti =  horzcat(ti, xi');
                        
                        % Hold on to the data lengths so we can keep track
                        % of which files are worst and should be ignored.
                        if isempty(find(isnan(zi(:,end)),1,'first'))
                            dataLengths = horzcat(dataLengths, size(zi(:,end),1));
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
                
                if size(dataStruct(velIndRef).z_corrected,2) > size(dataStruct(velIndRef).z_corrected,1)
                    dataStruct(length(Files)+j).z_average = dataStruct(velIndRef).z_corrected';
                else
                    dataStruct(length(Files)+j).z_average = dataStruct(velIndRef).z_corrected;
                end
                
                if size(dataStruct(velIndRef).d_corrected,2) > size(dataStruct(velIndRef).d_corrected,1)
                    dataStruct(length(Files)+j).d_average = dataStruct(velIndRef).d_corrected';
                else
                    dataStruct(length(Files)+j).d_average = dataStruct(velIndRef).d_corrected;
                end
                
                if size(dataStruct(velIndRef).t_approach,2) > size(dataStruct(velIndRef).t_approach,1)
                    dataStruct(length(Files)+j).t_average = dataStruct(velIndRef).t_approach';
                else
                    dataStruct(length(Files)+j).t_average = dataStruct(velIndRef).t_approach;
                end
                
                dataStruct(length(Files)+j).dt = dataStruct(velIndRef).dt;
                continue;
                
            end
            
            % Keep track of the shortest file
            [minLength,minInd] = min(dataLengths);
            [~,sortInd] = sort(dataLengths);
            [~,sortInd2] = sort(forceLimits);
            [~,sortInd3] = sort(timeLimits);
            dataStruct(length(Files)+j).shortestLength = minLength;
            dataStruct(length(Files)+j).shortestFile = currentFile(minInd);
            dataStruct(length(Files)+j).sortedFiles = currentFile(sortInd);
            dataStruct(length(Files)+j).sortedLengths = dataLengths(sortInd);
            dataStruct(length(Files)+j).sortedForceLimits = forceLimits(sortInd2);
            dataStruct(length(Files)+j).sortedForceFiles = currentFile(sortInd2);
            dataStruct(length(Files)+j).sortedTimeLimits = timeLimits(sortInd3);
            dataStruct(length(Files)+j).sortedTimeFiles = currentFile(sortInd3);
            
            % Average this load level's curves.
            % Interpolate according to the variance in the time array that
            % comes from averaging the results.
            t_interp = (1:size(ti,1))*median(dtVal);
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
                
            if size(z_interp,2) > size(z_interp,1)
                dataStruct(length(Files)+j).z_average = z_interp';
            else
                dataStruct(length(Files)+j).z_average = z_interp;
            end

            if size(d_interp,2) > size(d_interp,1)
                dataStruct(length(Files)+j).d_average = d_interp';
            else
                dataStruct(length(Files)+j).d_average = d_interp;
            end

            if size(t_interp,2) > size(t_interp,1)
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
    [z_max, z_max_ind] = max(dataStruct(length(Files)+i).z_average);
    
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
    if strcmp(filterType,'butter')
        % Create the butterworth
        [b,a] = butter(N,(cutoff_Hz)/(1/(2*dataStruct(length(Files)+i).dt)),'low'); % This makes a lowpass filter

        d_approach_smooth = (filter(b,a,dataStruct(length(Files)+i).d_approach)); % Next, apply the filter
    
    elseif strcmp(filterType,'IIR') || strcmp(filterType,'FIR')
        
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

    elseif strcmp(filterType,'movingAverageFFT')
        
        % Filtering using a moving average in the frequency domain
        f = linspace(0,1/(dataStruct(length(Files)+i).dt),length(dataStruct(length(Files)+i).d_approach));
        data = (fftshift(fft(dataStruct(length(Files)+i).d_approach)));
        M = movmean(data,7,'SamplePoints',f);
        
        figure
        ftBerk(dataStruct(length(Files)+i).t_approach,dataStruct(length(Files)+i).d_approach,1);
        hold on
        ftBerk(dataStruct(length(Files)+i).t_approach,abs(ifft(ifftshift(M))),1);
        hold off
        
        d_approach_smooth = abs(ifft(ifftshift(M)));
        delay = 0;
        
        % Force findrep to move in reverse...
        % Necessary or it won't find the right minimum
        findRep = 'reverse';
        
    elseif strcmp(filterType,'none')
        
        d_approach_smooth = dataStruct(length(Files)+i).d_approach;
        delay = 0;
        
    end
    
    if exist('sf','var')
        [~, dSmoothMin] = min(sf);
        d_approach_smooth = sf;
        z_approach_smooth = dataStruct(length(Files)+i).z_approach((delay+1):end);
    else
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
        tip_rep_pos = (length(tip_rep) - find(flipud(tip_rep)<0,1));       % Find last position above 0
        if isempty(tip_rep_pos) || tip_rep_pos == 0
            tip_rep_pos = 1;
        end
        
        tip_rep_pos_smooth = (length(tip_rep_smooth) - find(flipud(tip_rep_smooth)<0,1));       % Find last position above 0
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

%% Initialize
% Set the bounds of the for loop below so it either runs a single time (one
% instance of iterative or open-search fitting), or an array of times
% (multiple instances of iterative or open-search fitting).
if multiSim == 'y'
    nSims = n_terms_tot;
    voigtArray = 1:n_terms_tot;
else
    nSims = 1;
    voigtArray = n_terms_tot;
    skipNums = 0;
end

% Create timing Arrays
if timeOperation
    loopStartTimeStore = NaN(nSims,1);
    loopEndTimeStore = NaN(size(loopStartTimeStore));
    
    % Start the initial timer
    tic
end

%% Outer Loop
fclose('all');

% Make sure we are using the correct rows of the dataStruct with all of our
% experimental data. If we are using the 'a' (Averaged) analysis setting,
% we only care about the last few rows (one for each approach velocity used
% in the datasets). Alternately, if 'i' (Individual) is used, each row will
% need to be addressed SEPARATELY. So the indices will be rows 1 to
% n_datasets instead of n_datasets+j (where j is each separate approach
% velocity). For ease, the rows at the END of the dataStruct contain the
% averaged information.
switch avAnalysis
    case 'i'
        indShift = 0;
        loopInds = 1:length(Files);
    case 'a'
        indShift = length(Files);
        loopInds = 1:length(v_unique);
end

% Make all of our time storage arrays
if timeOperation
    fitSwitchStartTime = NaN(nSims,numel(loopInds));
    fitSwitchEndTime = NaN(size(fitSwitchStartTime));
    elasticParallelStartTime = NaN(nSims,numel(loopInds));
    elasticParallelEndTime = NaN(size(elasticParallelStartTime));
    elasticParallelBytes = NaN(size(elasticParallelEndTime));
    iterativeFitStartTime = NaN(nSims,numel(loopInds));
    iterativeFitEndTime = NaN(size(iterativeFitStartTime));
    parallelFitStartTime = NaN(nSims,numel(loopInds),nSims);
    parallelFitEndTime = NaN(size(parallelFitStartTime));
    parallelFitParallelBytes = NaN(size(parallelFitEndTime));
    parallelFitBreakdown = []; % Must be empty so it takes the correct datatype (Par object)
    parallelFitMemoryMonitor = cell(nSims,numel(loopInds));
end

%% Begin the Loop for our Model/Solver Configurations

if exist('dataPlot','var')
    clearvars dataPlot
end

for iSims = 1:nSims
    
    % Save the start time for this simulation loop
    if timeOperation
        loopStartTimeStore(iSims) = toc;
    end
    
    % Give the user feedback in the command window
    fprintf('\nRunning Simulation for %d Element(s)\n\n',voigtArray(iSims));
    
    % Grab the correct number of elements for this simulation loop
    n_terms = voigtArray(iSims);
    
    % Determine if the user asked us to skip this iteration.
    if ismember(n_terms,skipNums)
       fprintf('\nSkipped %d Retardance Terms per user request.\n\n',voigtArray(iSims));
       continue; 
    end

    % Open a Parallel Pool of MATLAB Workers
    if isempty(gcp('nocreate'))
        parpool('IdleTimeout', Inf)
    end

    % Start our text file log with the parameter information
    fid = fopen([pathname sprintf('/FitParams-nBranches_%d.txt',voigtArray(iSims))],'w');
    fprintf(fid,'Fitting Parameters, %s\r\n=================\r\n', date);
    fprintf(fid,'\r\nTip Radius: %4.3g\r\nVoigt Terms: %d\r\nScaling Factor: %4.3g',r_tip,n_terms,scaling);

    % Add the correct label to our text file, so someone looking later will
    % understand whether the files were analyzed separately or the average
    % was used.
    switch avAnalysis
        case 'i'
            fprintf(fid,'\r\nIndividual File Analysis\r\n=================\r\n');
        case 'a'
            fprintf(fid,'\r\nAveraged File Analysis\r\n=================\r\n');
    end
            
    % Loop through the distinct Approach Velocities (for avAnalysis=='a')
    % OR through every individual file (for avAnalysis=='i')
    for i_loop = loopInds
        
        % Grab the timestep and final time for use later
        dt = (dataStruct(indShift+i_loop).dt);
        st = dataStruct(indShift+i_loop).t_r(end);

        % Make the save label for our output files
        switch avAnalysis
            case 'i'
                saveLabel = sprintf('FileNo_%d',i_loop);
            case 'a'
                saveLabel = sprintf('LoadLevel%d',i_loop);
        end
        
        % Make the label for this iteration
        switch avAnalysis
            case 'i'
                fprintf(fid,'File No: %d\r\n',i_loop);
            case 'a'
                fprintf(fid,'Approach Velocity: %4.3g nm/s\r\n',v_unique(i_loop));
        end

        % Initialize the timing array for this loop
        if timeOperation
            timeStoreTemp = []; % Used to store the time for this loop
        end

        % Create our normalized data arrays (Lee and Radok). These are the
        % output quantities for each action integral, where the
        % normalization factor has been moved AWAY from the integral and
        % instead been applied to the output data (i.e. moved from the RHS
        % to the LHS of the equation).
        if smoothData == 'n'
            h_norm = ((8*sqrt(r_tip))/(3*(1-nu_sample))).*dataStruct(indShift+i_loop).h_r.^1.5;
            F_norm = ((3*(1-nu_sample))/(8*sqrt(r_tip))).*dataStruct(indShift+i_loop).F_r;
        else
            h_norm = ((8*sqrt(r_tip))/(3*(1-nu_sample))).*dataStruct(indShift+i_loop).h_r_smooth.^1.5;
            F_norm = ((3*(1-nu_sample))/(8*sqrt(r_tip))).*dataStruct(indShift+i_loop).F_r_smooth;
        end
        
        F_log = [];
        t_log = [];
        h_log = [];
        h_norm_log = [];

        % Get the input data organized
        if smoothData == 'n'
            F_log = log_scale(dataStruct(indShift+i_loop).F_r,dataStruct(indShift+i_loop).t_r,dt,st);
            t_log = log_scale(dataStruct(indShift+i_loop).t_r,dataStruct(indShift+i_loop).t_r,dt,st);
            h_log = log_scale(dataStruct(indShift+i_loop).h_r,dataStruct(indShift+i_loop).t_r,dt,st);

            negativeIndexes = dataStruct(indShift+i_loop).F_r <=  0;
            x_fit = [dataStruct(indShift+i_loop).t_r, dataStruct(indShift+i_loop).F_r];
            x_fit(negativeIndexes,2) = 1e-25;
            timeVec = dataStruct(indShift+i_loop).t_r;
            maxwellInputs = [dataStruct(indShift+i_loop).t_r(dataStruct(indShift+i_loop).t_r>=0), (dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>=0).^1.5)];
            maxwellTimeVec = dataStruct(indShift+i_loop).t_r(dataStruct(indShift+i_loop).t_r>=0);
        else
            F_log = log_scale(dataStruct(indShift+i_loop).F_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);
            t_log = log_scale(dataStruct(indShift+i_loop).t_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);
            h_log = log_scale(dataStruct(indShift+i_loop).h_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);

            negativeIndexes = dataStruct(indShift+i_loop).F_r_smooth <=  0;
            x_fit = [dataStruct(indShift+i_loop).t_r_smooth, dataStruct(indShift+i_loop).F_r_smooth];
            x_fit(negativeIndexes,2) = 1e-25;
            timeVec = dataStruct(indShift+i_loop).t_r_smooth;
            maxwellInputs = [dataStruct(indShift+i_loop).t_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0), (dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0).^1.5)];
            maxwellTimeVec = dataStruct(indShift+i_loop).t_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0);
        end

        dataStruct(indShift+i_loop).h_norm = h_norm;
        dataStruct(indShift+i_loop).F_norm = F_norm;

        dataStruct(indShift+i_loop).h_log = h_log;
        h_norm_log = log_scale(h_norm,timeVec,dt,st);
        dataStruct(indShift+i_loop).h_norm_log = h_norm_log;
        F_norm_log = log_scale(F_norm,timeVec,dt,st);
        dataStruct(indShift+i_loop).F_norm_log = F_norm_log;
        
        % Open Parallel Pool of MATLAB Workers
        if isempty(gcp('nocreate'))
            parpool('IdleTimeout', Inf)
        end

        if timeOperation
            fitSwitchStartTime(iSims,i_loop) = toc;
        end

        % Jump into the switch-case that controls whether the "open" or
        % "iterative" parameter search method is used.
        switch lower(fitMethod)
            
            case 'iterative'
                % Initialize the cell arrays that will store the parameters
                % from every iteration up to i_loop
                maxwellParams_loop = {};
                convParams_loop = {};

                if strcmp(elasticSetting,'y')
                    % Get the best elastic parameter for our iterative fit
                    [~,~,~,lb_elastic,ub_elastic,~,~] = makeGeneralizedVoigtModel('y','n',0,timeVec,dt,st,minTimescale,forwardFitTimescale);

                    x_fit_elastic = [log_scale(x_fit(:,1),x_fit(:,1),dt,st); log_scale(x_fit(:,2),x_fit(:,1),dt,st)]';
                    F_conv_wrapper_elastic = @(c,inputs) c(1).*inputs(:,2);

                    y_fit_elastic = (h_norm_log');

                    if fitElasticFirst
                        beta_dist_elastic = zeros(length(ub_elastic),n_samples);
                        beta0_dist_elastic = zeros(size(beta_dist_elastic));
                        resnorm_dist_elastic = zeros(1,n_samples);

                        if timeOperation
                            elasticParallelStartTime(iSims,i_loop) = toc;
                            ticBytes(gcp);
                        end

                        ub_rand_elastic = log10(ub_elastic(1))-1;
                        lb_rand_elastic = log10(lb_elastic(1))+1;
                        beta0_elastic_array = logspace(ub_rand_elastic,lb_rand_elastic,n_samples);

                        progressString = sprintf('Investigating Elastic Parameter\nParallel Search Running...');
                        hbar = parfor_progressbar(n_samples,progressString);
                        warning('off');
                        
                        parfor i = 1:n_samples
                            lsqoptions = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                                'MaxFunctionEvaluations', n_maxIterations,...
                                'MaxIterations', n_maxIterations,...
                                'FiniteDifferenceType','central',...
                                'FunctionTolerance', 0,...
                                'OptimalityTolerance', 0,...
                                'StepTolerance', scaling,...
                                'Display', 'none');

                            [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                lsqcurvefit(F_conv_wrapper_elastic,beta0_elastic_array(i),x_fit_elastic,y_fit_elastic,lb_elastic,ub_elastic,lsqoptions);

                            beta_dist_elastic(:,i) = betatemp;
                            beta0_dist_elastic(:,i) = beta0_elastic_array(i);
                            resnorm_dist_elastic(i) = sum(((F_conv_wrapper_elastic(betatemp,x_fit_elastic)-y_fit_elastic).^2)./movvar(y_fit_elastic,3))./(length(y_fit_elastic)-length(betatemp)); % Mean Squared Error (MSE)
                            
                            hbar.iterate(1) % Increase progressbar by 1 iteration
                        end

                        if timeOperation
                            elasticParallelEndTime(iSims,i_loop) = toc;
                            elasticParallelBytes(iSims,i_loop) = getfield(tocBytes(gcp),{1});

                            fprintf('\nThe Elastic Parameter parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                                (elasticParallelEndTime(iSims,i_loop)-elasticParallelStartTime(iSims,i_loop))/60, elasticParallelBytes(iSims,i_loop)*1e-6);
                        end
                        
                        close(hbar)
                        warning('on');

                        resnorm_dist_elastic_temp = resnorm_dist_elastic;
                        resnorm_dist_elastic_temp(resnorm_dist_elastic_temp == 0) = NaN;
                        [~,idx] = min(resnorm_dist_elastic_temp,[],'omitnan');
                        bestElasticTerm = beta_dist_elastic(:,idx);

                    else
                        if timeOperation
                            elasticParallelStartTime(iSims,i_loop) = toc;
                            ticBytes(gcp);
                            elasticParallelEndTime(iSims,i_loop) = toc;
                            elasticParallelBytes(iSims,i_loop) = getfield(tocBytes(gcp),{1});
                        end
                    end
                end

                if timeOperation
                    iterativeFitStartTime(iSims,i_loop) = toc;
                end

                % Perform the iterative fitting process for this loop
                for i = 1:n_terms

                    % We need to grab the correct size for our parameter
                    % set in this iteration. We will create the model
                    % functions below.
                    [~,~,~,~,ub_loop,~,~] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,i,timeVec,dt,st,minTimescale,forwardFitTimescale);
                    [~,~,~,~,padSize_loop] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellTimeVec,dt,st);
                    [~,~,~,~,~,ub_loop_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellTimeVec,minTimescale,padSize_loop,forwardFitTimescale);

                    % Get the indices of the Moduli and Characteristic
                    % Times for each model.
                    [tauInds,modulusInds] = getParamIndices(ub_loop,elasticSetting,fluidSetting);
                    [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(ub_loop_maxwell,elasticSetting_maxwell,fluidSetting_maxwell);

                    % Add the elastic term to our list, if used
                    if strcmp(elasticSetting,'y')
                        elasticInd = 1;
                        modulusInds = horzcat(elasticInd,modulusInds);
                    else
                        elasticInd = 0;
                    end

                    % Add the elastic term to our list, if used
                    if strcmp(elasticSetting_maxwell,'y')
                        elasticInd = 1;
                        modulusInds_maxwell = horzcat(elasticInd,modulusInds_maxwell);
                    else
                        elasticInd = 0;
                    end

                    % Get our observables for this loop. If previously
                    % enabled, this will clip the datasets to a reasonable
                    % bound for the current model configuration.
                    if strcmp(smoothData,'n')
                        if clipData
                            if fitLog
                                y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds))),timeVec,dt,st);
                                y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec,dt,st);
                            else
                                y_fit_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                            end

                            maxwellInputs_loop = maxwellInputs(x_fit(:,1) < max(ub_loop_maxwell(tauInds)),:);
                            x_fit_loop = x_fit(x_fit(:,1) < max(ub_loop(tauInds)),:);
                            t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                            h_r_loop = dataStruct(indShift+i_loop).h_r(x_fit(:,1) < max(ub_loop(tauInds)));
                            h_norm_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                            F_r_loop = dataStruct(indShift+i_loop).F_r(x_fit(:,1) < max(ub_loop(tauInds)));
                            F_norm_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                            h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                            F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                        else
                            if fitLog
                                y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm,timeVec,dt,st);
                                y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm,timeVec,dt,st);
                            else
                                y_fit_loop = dataStruct(indShift+i_loop).h_norm;
                                y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm;
                            end

                            maxwellInputs_loop = maxwellInputs;
                            x_fit_loop = x_fit;
                            t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                            h_r_loop = dataStruct(indShift+i_loop).h_r;
                            h_norm_loop = dataStruct(indShift+i_loop).h_norm;
                            F_r_loop = dataStruct(indShift+i_loop).F_r;
                            F_norm_loop = dataStruct(indShift+i_loop).F_norm;
                            h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                            F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                        end
                    else
                        if clipData
                            if fitLog
                                y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds))),timeVec,dt,st);
                                y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec,dt,st);
                            else
                                y_fit_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                            end

                            maxwellInputs_loop = maxwellInputs(x_fit(:,1) < max(ub_loop_maxwell(tauInds)),:);
                            x_fit_loop = x_fit(x_fit(:,1) < max(ub_loop(tauInds)),:);
                            t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                            h_r_loop = dataStruct(indShift+i_loop).h_r_smooth(x_fit(:,1) < max(ub_loop(tauInds)));
                            h_norm_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                            F_r_loop = dataStruct(indShift+i_loop).F_r_smooth(x_fit(:,1) < max(ub_loop(tauInds)));
                            F_norm_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                            h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                            F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                        else
                            if fitLog
                                y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm,timeVec,dt,st);
                                y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm,timeVec,dt,st);
                            else
                                y_fit_loop = dataStruct(indShift+i_loop).h_norm;
                                y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm;
                            end

                            maxwellInputs_loop = maxwellInputs;
                            x_fit_loop = x_fit;
                            t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                            h_r_loop = dataStruct(indShift+i_loop).h_r_smooth;
                            h_norm_loop = dataStruct(indShift+i_loop).h_norm;
                            F_r_loop = dataStruct(indShift+i_loop).F_r_smooth;
                            F_norm_loop = dataStruct(indShift+i_loop).F_norm;
                            h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                            F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                        end
                    end

                    % Now we create our function definitions for fitting.
                    % Make the Generalized Voigt Model for this loop
                    [~,F_conv_loop,F_conv_wrapper_loop,lb_loop,ub_loop,~,~] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,i,x_fit_loop(:,1),dt,st,minTimescale,forwardFitTimescale);

                    % Create the Gen. Maxwell Force Convolution Model
                    [~,~,h_conv_maxwell_loop,h_conv_maxwell_wrapper_loop,padSize_loop] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellInputs_loop(:,1),dt,st);

                    % Get the upper and lower bounds for the maxwell model
                    [~,~,~,~,lb_loop_maxwell,ub_loop_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellInputs_loop(:,1),minTimescale,padSize_loop,forwardFitTimescale);

                    % Initialize our current parameter set. Values that
                    % remain NaN when these arrays get to the parfor loop
                    % are REPLACED inside the loop. We use NaN to identify
                    % the newly-introduced model parameters so they can be
                    % randomized/updated inside of the parfor, and we keep
                    % the old parameters (from the last loop) intact.
                    beta0_maxwell = NaN(size(ub_loop_maxwell));
                    beta0_conv = NaN(size(ub_loop));

                    % These are the "reasonable" bounds for our guesses.
                    % They are given in log10 form (i.e. -1 is 10^(-1) when
                    % generating the random parameters).
                    ub_randLimit = 1*ones(size(ub_loop));
                    ub_randLimit(tauInds) = ub_loop(tauInds);
                    lb_randLimit = -5*ones(size(lb_loop));
                    lb_randLimit(tauInds) = lb_loop(tauInds);

                    ub_randLimit_maxwell = 5*ones(size(ub_loop_maxwell));
                    ub_randLimit_maxwell(tauInds_maxwell) = ub_loop_maxwell(tauInds_maxwell);
                    lb_randLimit_maxwell = -1*ones(size(lb_loop_maxwell));
                    lb_randLimit_maxwell(tauInds_maxwell) = lb_loop_maxwell(tauInds_maxwell);

                    % For the Voigt model, we need to bring our previous
                    % iteration parameters into the initial guess array we
                    % initialized above.
                    if i > 1
                        if fluidSetting == 'y'
                            beta0_conv(1:length(convParams_loop{i-1})-1) = convParams_loop{i-1}(1:end-1);
                            beta0_conv(end) = 0;
                        else
                            beta0_conv(1:length(convParams_loop{i-1})) = convParams_loop{i-1};
                            if strcmp(elasticSetting,'y') && isnan(beta0_conv(1))
                                beta0_conv(1) = bestElasticTerm;
                            end
                        end
                    else
                        % Manual override beta0 for random sampling
                        beta0_conv = ((ub_loop-lb_loop).*rand(size(ub_loop)) + lb_loop);
                    end

                    % For the Maxwell model, we need to bring our previous
                    % iteration parameters into the initial guess array we
                    % initialized above, as well.
                    if i > 1
                        if fluidSetting_maxwell == 'y'
                            beta0_maxwell(1:length(maxwellParams_loop{i-1})-1) = maxwellParams_loop{i-1}(1:end-1);
                            beta0_maxwell(end) = 0;
                        else
                            beta0_maxwell(1:length(maxwellParams_loop{i-1})) = maxwellParams_loop{i-1};
                            if strcmp(elasticSetting_maxwell,'y') && isnan(beta0_maxwell(1))
                                beta0_maxwell(1) = 1/bestElasticTerm;
                            end
                        end
                    else
                        % Manual override beta0 for random sampling
                        beta0_maxwell = ((ub_loop_maxwell-lb_loop_maxwell).*rand(size(ub_loop_maxwell)) + lb_loop_maxwell);
                    end
                    
                    % Make our storage arrays. These MUST be preallocated
                    % or the parfor loop will complain.
                    beta_dist = zeros(length(ub_loop),n_samples);
                    beta0_dist = zeros(size(beta_dist));
                    resnorm_dist = zeros(1,n_samples);
                    beta_dist_maxwell = zeros(length(ub_loop_maxwell),n_samples);
                    beta0_dist_maxwell = zeros(size(beta_dist_maxwell));
                    resnorm_dist_maxwell = zeros(1,n_samples);

                    % Start the timer and make our timing placeholders.
                    % Again, these must be preallocated or the parfor loop
                    % will break.
                    if timeOperation
                        parallelFitStartTime(iSims,i_loop,i) = toc;
                        convTime = Par(n_samples);
                        maxwellTime = Par(n_samples);
                        if i == 1
                           convTimeTotal = [];
                           maxwellTimeTotal = [];
                        end
                        ticBytes(gcp);
                        if ~exist('userMem','var')
                            [userMem,systemMem] = memory;
                        end
                        ramused = NaN(1,n_samples);
                        ramav = NaN(1,n_samples);
                        ramtime = NaN(1,n_samples);
                    end
                    
                    progressString = sprintf('%d-Term Viscoelastic Models, Iterative Search\nInvestigating Material Parameters\nParallel Search Running...',i);
                    hbar = parfor_progressbar(n_samples,progressString);
                    warning('off');
                    
                    % Parallel Loop for Random Search
                    parfor j = 1:n_samples
                        warning('off');
                        lastwarn(''); % Clear last warning

                        try
                            % Initialize
                            newInds_conv = [];
                            newInds_maxwell = [];
                            beta0_conv_loop = [];
                            beta0_maxwell_loop = [];

                            % Find the new indices
                            newInds_conv = isnan(beta0_conv);
                            newInds_maxwell = isnan(beta0_maxwell);

                            % Grab the base values for the beta0 cases
                            beta0_conv_loop = beta0_conv;
                            beta0_maxwell_loop = beta0_maxwell;

                            % Control the bounds AND the initial parameters
                            ub_fluidityRandLimit_init = 1e-20;

                            % Include or remove elastic term from bound
                            % loosening that happens below. Use 1 for
                            % inclusion, empty brackets [] for ignoring.
                            elasticInd = 1;

                            % Pick the right functions to use
                            if fitLog
                                F_conv_func = @(c,inputs) F_conv_wrapper_loop(c,inputs);
                                h_maxwell_func = @(c,inputs) h_conv_maxwell_wrapper_loop(c,inputs);
                            else
                                F_conv_func = @(c,inputs) F_conv_loop(c,inputs);
                                h_maxwell_func = @(c,inputs) h_conv_maxwell_loop(c,inputs);
                            end

                            % Settings for lsqcurvefit
                            lsqoptions = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                                'MaxFunctionEvaluations', n_maxIterations,...
                                'MaxIterations', n_maxIterations,...
                                'FiniteDifferenceType','central',...
                                'FunctionTolerance', scaling,...
                                'OptimalityTolerance', 0,...
                                'StepTolerance', scaling,...
                                'Display', 'none');

                            % Create random starting point
                            [beta0_conv_loop,~,~] = makeRandomParams(beta0_conv,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,newInds_conv);

                            % Create random starting point
                            [beta0_maxwell_loop,~,~] = makeRandomParams(beta0_maxwell,ub_randLimit_maxwell,lb_randLimit_maxwell,elasticSetting_maxwell,fluidSetting_maxwell,newInds_maxwell);

                            % One last check of the parameters to try...
                            if any(isnan(beta0_conv_loop))
                                tempInd = isnan(beta0_conv_loop);
                                beta0_conv_loop(tempInd) = lb_loop(tempInd)+rand.*(ub_loop(tempInd)-lb_loop(tempInd));
                            end

                            if any(isnan(beta0_maxwell_loop))
                                tempInd = isnan(beta0_maxwell_loop);
                                beta0_maxwell_loop(tempInd) = lb_loop_maxwell(tempInd)+rand.*(ub_loop_maxwell(tempInd)-lb_loop_maxwell(tempInd));
                            end

                            if timeOperation
                                Par.tic;
                            end

                            [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                lsqcurvefit(F_conv_func,beta0_conv_loop,x_fit_loop,y_fit_loop,lb_loop,ub_loop,lsqoptions);

                            if timeOperation
                                convTime(j) = Par.toc;
                            end

                            beta_dist(:,j) = betatemp';
                            beta0_dist(:,j) = beta0_conv_loop;
                            resnorm_dist(j) = sum(((F_conv_func(betatemp,x_fit_loop)-y_fit_loop).^2)./movvar(y_fit_loop,3))./(length(y_fit_loop)-length(betatemp)); % Mean Squared Error (MSE)

                            if timeOperation
                                Par.tic;
                            end

                            [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                lsqcurvefit(h_maxwell_func,beta0_maxwell_loop,maxwellInputs_loop,y_fit_maxwell_loop,lb_loop_maxwell,ub_loop_maxwell,lsqoptions);

                            if timeOperation
                                maxwellTime(j) = Par.toc;
                            end

                            beta_dist_maxwell(:,j) = betatemp';
                            beta0_dist_maxwell(:,j) = beta0_maxwell_loop;
                            resnorm_dist_maxwell(j) = sum(((h_maxwell_func(betatemp,maxwellInputs_loop)-y_fit_maxwell_loop).^2)./movvar(y_fit_maxwell_loop,3))./(length(y_fit_maxwell_loop)-length(betatemp)); % Mean Squared Error (MSE)

                            if timeOperation
                                % Track memory
                                ramused(j) = userMem.MemUsedMATLAB;
                                ramav(j) = systemMem.PhysicalMemory.Available;
                                ramtime(j) = now;
                            end

                        catch ERROR

                            beta_dist(:,j) = NaN;
                            beta0_dist(:,j) = NaN;
                            resnorm_dist(j) = NaN;

                            beta_dist_maxwell(:,j) = NaN;
                            beta0_dist_maxwell(:,j) = NaN;
                            resnorm_dist_maxwell(j) = NaN;

                            if timeOperation
                                convTime(j) = 1;
                                maxwellTime(j) = 1;
                            end

                            voigtData = [0,0];
                            maxwellData = [0,0];

                            fprintf('\nTerm Number %d Search: Entered a NaN Row for loop #%d because there was an ERROR inside the parfor loop:\n',i,j);
                            disp(ERROR.message)
                        end
                        hbar.iterate(1) % Increase progressbar by 1 iteration
                    end

                    if timeOperation
                        parallelFitEndTime(iSims,i_loop,i) = toc;
                        parallelFitParallelBytes(iSims,i_loop,i) = getfield(tocBytes(gcp),{1});

                        stop(convTime);
                        stop(maxwellTime);

                        convTimeTotal = [convTimeTotal convTime];
                        maxwellTimeTotal = [maxwellTimeTotal maxwellTime];

                        fprintf('\nThe parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                            (parallelFitEndTime(iSims,i_loop,i)-parallelFitStartTime(iSims,i_loop,i))/60, parallelFitParallelBytes(iSims,i_loop,i)*1e-6);
                    end
                    
                    close(hbar)
                    warning('on');

                    resnorm_dist_temp = resnorm_dist;
                    resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                    resnorm_dist_temp( any(abs(log10(beta_dist(modulusInds,:)) - log10(repmat(lb_loop(modulusInds)',1,size(beta_dist(modulusInds,:),2))))<tossThresh,1) ) = NaN;
                    resnorm_dist_temp( any(abs(log10(repmat(ub_loop(modulusInds)',1,size(beta_dist(modulusInds,:),2))) - log10(beta_dist(modulusInds,:)))<tossThresh,1) ) = NaN;
                    [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                    convParams_loop{i} = beta_dist(:,idx);

                    resnorm_dist_maxwell_temp = resnorm_dist_maxwell;
                    resnorm_dist_maxwell_temp(resnorm_dist_maxwell_temp == 0) = NaN;
                    resnorm_dist_maxwell_temp( any(abs(log10(beta_dist_maxwell(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_maxwell(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                    resnorm_dist_maxwell_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_maxwell(modulusInds_maxwell,:),2))) - log10(beta_dist_maxwell(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                    [~,idx_maxwell] = min(resnorm_dist_maxwell_temp);
                    maxwellParams_loop{i} = beta_dist_maxwell(:,idx_maxwell);

                end

                if timeOperation
                    iterativeFitEndTime(iSims,i_loop) = toc;
                    parallelFitBreakdown{iSims,i_loop} = {convTimeTotal,maxwellTimeTotal};
                    parallelFitMemoryMonitor{iSims,i_loop} = vertcat(ramused,ramav,ramtime);
                end

                dataStruct(indShift+i_loop).convParams_loop = convParams_loop;
                dataStruct(indShift+i_loop).maxwellParams_loop = maxwellParams_loop;

                n_ranks = 10; % Number of places to save the parameter sets of

                resnorm_dist_temp = resnorm_dist;
                resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                resnorm_dist_temp( any(abs(log10(beta_dist(modulusInds,:)) - log10(repmat(lb_loop(modulusInds)',1,size(beta_dist(modulusInds,:),2))))<tossThresh,1) ) = NaN;
                resnorm_dist_temp( any(abs(log10(repmat(ub_loop(modulusInds)',1,size(beta_dist(modulusInds,:),2))) - log10(beta_dist(modulusInds,:)))<tossThresh,1) ) = NaN;
                [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                dataStruct(indShift+i_loop).convParams = beta_dist(:,idx);

                [~,I] = sort(resnorm_dist_temp,'MissingPlacement','last');
                dataStruct(indShift+i_loop).convParams_alternatives = beta_dist(:,I(1:n_ranks));

                resnorm_dist_maxwell_temp = resnorm_dist_maxwell;
                resnorm_dist_maxwell_temp(resnorm_dist_maxwell_temp == 0) = NaN;
                resnorm_dist_maxwell_temp( any(abs(log10(beta_dist_maxwell(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_maxwell(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                resnorm_dist_maxwell_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_maxwell(modulusInds_maxwell,:),2))) - log10(beta_dist_maxwell(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                [~,idx_maxwell] = min(resnorm_dist_maxwell_temp);
                dataStruct(indShift+i_loop).maxwellParams = beta_dist_maxwell(:,idx_maxwell);

                [~,I_maxwell] = sort(resnorm_dist_maxwell_temp,'MissingPlacement','last');
                dataStruct(indShift+i_loop).maxwellParams_alternatives = beta_dist_maxwell(:,I_maxwell(1:n_ranks));

                % Print Parameters
                fprintf(fid,'\r\n\r\nGeneralized Kelvin-Voigt\r\n');
                if elasticSetting == 'y'
                    fprintf(fid,'J_g = %3.6g [1/Pa]\r\n',beta_dist(1,idx));

                    for iii = 1:n_terms
                        fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2,idx),iii,beta_dist(iii*2+1,idx));
                    end
                    if fluidSetting == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,idx));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                else
                    fprintf(fid,'No elastic response used.\r\n');

                    for iii = 1:n_terms
                        fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2-1,idx),iii,beta_dist(iii*2,idx));
                    end
                    if fluidSetting == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,idx));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                end

                fprintf(fid,'\r\n\r\nGeneralized Maxwell\r\n');
                if elasticSetting_maxwell == 'y'
                    fprintf(fid,'E_e = %3.6g [Pa]\r\n',beta_dist_maxwell(1,idx_maxwell));

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,idx_maxwell),iii,beta_dist_maxwell(iii*2+1,idx_maxwell));
                    end
                    if fluidSetting_maxwell == 'y'
                        fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                else
                    fprintf(fid,'No elastic response used.\r\n');

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,idx_maxwell),iii,beta_dist_maxwell(iii*2,idx_maxwell));
                    end
                    if fluidSetting_maxwell == 'y'
                        fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                end

                % Print alternate parameters
                fprintf(fid,'\r\n\r\nAlternative Fitting Parameters (2:10)\r\n\r\n');
                for alti = 2:n_ranks
                    fprintf(fid,' \r\n\r\nRank Number %d:',alti);

                    % Print GKV
                    fprintf(fid,'\r\n\r\nGeneralized Kelvin-Voigt\r\n');
                    if elasticSetting == 'y'
                        fprintf(fid,'J_g = %3.6g [1/Pa]\r\n',beta_dist(1,I(alti)));

                        for iii = 1:n_terms
                            fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2,I(alti)),iii,beta_dist(iii*2+1,I(alti)));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,I(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2-1,I(alti)),iii,beta_dist(iii*2,I(alti)));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,I(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                    % Print Maxwell
                    fprintf(fid,'\r\n\r\nGeneralized Maxwell\r\n');
                    if elasticSetting_maxwell == 'y'
                        fprintf(fid,'E_e =%3.6g [Pa]\r\n',beta_dist_maxwell(1,I_maxwell(alti)));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2+1,I_maxwell(alti)));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2,I_maxwell(alti)));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                end

            case 'open'

                if timeOperation
                    iterativeFitStartTime(iSims,i_loop) = toc;
                end

                % Make the Generalized Voigt for this loop
                [U_func,F_conv,F_conv_wrapper,lb,ub,subref,selector] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,timeVec,dt,st,minTimescale,forwardFitTimescale);

                % Create the Gen. Maxwell Force Convolution Model
                [G_func_maxwell,Q_func_maxwell,h_conv_maxwell,h_conv_maxwell_wrapper,padSize] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_terms,maxwellTimeVec,dt,st);

                % Create the Harmonic Gen. Maxwell Models
                [E_storage,E_loss,lossAngle,harmonic_wrapper,lb_maxwell,ub_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_terms,maxwellTimeVec,minTimescale,padSize,forwardFitTimescale);

                [tauInds,modulusInds] = getParamIndices(ub,elasticSetting,fluidSetting);
                [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(ub_maxwell,elasticSetting_maxwell,fluidSetting_maxwell);

                if strcmp(elasticSetting,'y')
                    elasticInd = 1;
                    modulusInds = horzcat(elasticInd,modulusInds);
                else
                    elasticInd = 0;
                end

                if strcmp(elasticSetting_maxwell,'y')
                    elasticInd = 1;
                    modulusInds_maxwell = horzcat(elasticInd,modulusInds_maxwell);
                else
                    elasticInd = 0;
                end
                
                % These are the "reasonable" bounds for our guesses.
                % They are given in log10 form (i.e. -1 is 10^(-1) when
                % generating the random parameters).
                ub_randLimit = 1*ones(size(ub));
                ub_randLimit(tauInds) = ub(tauInds);
                lb_randLimit = -5*ones(size(lb));
                lb_randLimit(tauInds) = lb(tauInds);

                ub_randLimit_maxwell = 5*ones(size(ub_maxwell));
                ub_randLimit_maxwell(tauInds_maxwell) = ub_maxwell(tauInds_maxwell);
                lb_randLimit_maxwell = -1*ones(size(lb_maxwell));
                lb_randLimit_maxwell(tauInds_maxwell) = lb_maxwell(tauInds_maxwell);

                beta_dist = zeros(length(ub),n_samples);
                beta0_dist = zeros(size(beta_dist));
                resnorm_dist = zeros(1,n_samples);
                beta_dist_maxwell = zeros(length(ub_maxwell),n_samples);
                beta0_dist_maxwell = zeros(size(beta_dist_maxwell));
                resnorm_dist_maxwell = zeros(1,n_samples);
                
                if fitLog
                    y_fit = dataStruct(indShift+i_loop).h_norm_log;
                    y_fit_maxwell = dataStruct(indShift+i_loop).F_norm_log;
                else
                    y_fit = dataStruct(indShift+i_loop).h_norm;
                    y_fit_maxwell = dataStruct(indShift+i_loop).F_norm;
                end

                if timeOperation
                    parallelFitStartTime(iSims,i_loop,n_terms) = toc;
                    convTime = Par(n_samples);
                    maxwellTime = Par(n_samples);
                    if i == 1
                       convTimeTotal = [];
                       maxwellTimeTotal = [];
                    end
                    ticBytes(gcp);
                    if ~exist('userMem','var')
                        [userMem,systemMem] = memory;
                    end
                    ramused = NaN(1,n_samples);
                    ramav = NaN(1,n_samples);
                    ramtime = NaN(1,n_samples);
                end
                
                progressString = sprintf('%d-Term Viscoelastic Models, Open Search\nInvestigating Material Parameters\nParallel Search Running...',n_terms);
                hbar = parfor_progressbar(n_samples,progressString);
                warning('off');

                % Parallel Loop for Random Search
                parfor j = 1:n_samples
                    warning('off');
                    lastwarn(''); % Clear last warning

                    try
                        % Initialize
                        beta0_conv = [];
                        beta0_maxwell = [];

                        beta0_maxwell = NaN(size(ub_maxwell));
                        beta0_conv = NaN(size(ub));

                        ub_fluidityRandLimit_init = 1e-20;

                        % Settings for lsqcurvefit
                        lsqoptions = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                            'MaxFunctionEvaluations', n_maxIterations,...
                            'MaxIterations', n_maxIterations,...
                            'FiniteDifferenceType','central',...
                            'FunctionTolerance', scaling,...
                            'OptimalityTolerance', scaling,...
                            'StepTolerance', scaling,...
                            'Display', 'none');

                        % Pick the right functions to use
                        if fitLog
                            F_conv_func = @(c,inputs) F_conv_wrapper(c,inputs);
                            h_maxwell_func = @(c,inputs) h_conv_maxwell_wrapper(c,inputs);
                        else
                            F_conv_func = @(c,inputs) F_conv(c,inputs);
                            h_maxwell_func = @(c,inputs) h_conv_maxwell(c,inputs);
                        end

                        % Create random starting point
                        [beta0_conv_loop,~,~] = makeRandomParams(beta0_conv,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,isnan(beta0_conv));

                        % Create random starting point
                        [beta0_maxwell_loop,~,~] = makeRandomParams(beta0_maxwell,ub_randLimit_maxwell,lb_randLimit_maxwell,elasticSetting_maxwell,fluidSetting_maxwell,isnan(beta0_maxwell));

                        if timeOperation
                            Par.tic;
                        end

                        [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                            lsqcurvefit(F_conv_func,beta0_conv_loop,x_fit,y_fit,lb,ub,lsqoptions);

                        if timeOperation
                            convTime(j) = Par.toc;
                        end

                        beta_dist(:,j) = betatemp';
                        beta0_dist(:,j) = beta0_conv_loop;
                        resnorm_dist(j) = sum(((F_conv_func(betatemp,x_fit)-y_fit).^2)./movvar(y_fit,3))./(length(y_fit)-length(betatemp)); % Mean Squared Error (MSE)

                        if timeOperation
                            Par.tic;
                        end

                        [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                    lsqcurvefit(h_maxwell_func,beta0_maxwell_loop,maxwellInputs,y_fit_maxwell,lb_maxwell,ub_maxwell,lsqoptions);

                        if timeOperation
                            maxwellTime(j) = Par.toc;
                        end

                        beta_dist_maxwell(:,j) = betatemp';
                        beta0_dist_maxwell(:,j) = beta0_maxwell_loop;
                        resnorm_dist_maxwell(j) = sum(((h_maxwell_func(betatemp,maxwellInputs)-y_fit_maxwell).^2)./movvar(y_fit_maxwell,3))./(length(y_fit_maxwell)-length(betatemp)); % Mean Squared Error (MSE)

                        if timeOperation
                            % Track memory
                            ramused(j) = userMem.MemUsedMATLAB;
                            ramav(j) = systemMem.PhysicalMemory.Available;
                            ramtime(j) = now;
                        end

                    catch ERROR
                        beta_dist(:,j) = NaN;
                        beta0_dist(:,j) = NaN;
                        resnorm_dist(j) = NaN;

                        beta_dist_maxwell(:,j) = NaN;
                        beta0_dist_maxwell(:,j) = NaN;
                        resnorm_dist_maxwell(j) = NaN;

                        if timeOperation
                            convTime(j) = 1;
                            maxwellTime(j) = 1;
                        end

                        voigtData = [0,0];
                        maxwellData = [0,0];

                        fprintf('\nTerm Number %d Search: Entered a NaN Row for loop #%d because there was an ERROR inside the parfor loop:\n',n_terms,j);
                        disp(ERROR.message)
                    end
                    hbar.iterate(1) % Increase progressbar by 1 iteration
                end

                if timeOperation
                    parallelFitEndTime(iSims,i_loop,n_terms) = toc;
                    parallelFitParallelBytes(iSims,i_loop,n_terms) = getfield(tocBytes(gcp),{1});

                    stop(convTime);
                    stop(maxwellTime);

                    convTimeTotal = convTime;
                    maxwellTimeTotal = maxwellTime;

                    fprintf('\nThe parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                        (parallelFitEndTime(iSims,i_loop,n_terms)-parallelFitStartTime(iSims,i_loop,n_terms))/60, parallelFitParallelBytes(iSims,i_loop,n_terms)*1e-6);
                end
                
                close(hbar)
                warning('on');

                n_ranks = 10; % Number of places to save the parameter sets of

                resnorm_dist_temp = resnorm_dist;
                resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                dataStruct(indShift+i_loop).convParams = beta_dist(:,idx);

                [~,I] = sort(resnorm_dist_temp,'MissingPlacement','last');
                dataStruct(indShift+i_loop).convParams_alternatives = beta_dist(:,I(1:n_ranks));

                resnorm_dist_maxwell_temp = resnorm_dist_maxwell;
                resnorm_dist_maxwell_temp(resnorm_dist_maxwell_temp == 0) = NaN;
                [~,idx_maxwell] = min(resnorm_dist_maxwell_temp);
                dataStruct(indShift+i_loop).maxwellParams = beta_dist_maxwell(:,idx_maxwell);

                [~,I_maxwell] = sort(resnorm_dist_maxwell_temp,'MissingPlacement','last');
                dataStruct(indShift+i_loop).maxwellParams_alternatives = beta_dist_maxwell(:,I_maxwell(1:n_ranks));

                if timeOperation
                    iterativeFitEndTime(iSims,i_loop) = toc;
                    parallelFitBreakdown{iSims,i_loop} = {convTimeTotal,maxwellTimeTotal,schaperyTimeTotal};
                    parallelFitMemoryMonitor{iSims,i_loop} = vertcat(ramused,ramav,ramtime);
                end

                % Print Parameters to our output file in a nicer, formatted
                % way.
                fprintf(fid,'\r\n\r\nGeneralized Kelvin-Voigt\r\n');
                if elasticSetting == 'y'
                    fprintf(fid,'J_g = %3.6g [1/Pa]\r\n',beta_dist(1,idx));

                    for iii = 1:n_terms
                        fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2,idx),iii,beta_dist(iii*2+1,idx));
                    end
                    if fluidSetting == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,idx));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                else
                    fprintf(fid,'No elastic response used.\r\n');

                    for iii = 1:n_terms
                        fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2-1,idx),iii,beta_dist(iii*2,idx));
                    end
                    if fluidSetting == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,idx));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                end

                fprintf(fid,'\r\n\r\nGeneralized Maxwell\r\n');
                if elasticSetting == 'y'
                    fprintf(fid,'E_e = %3.6g [Pa]\r\n',beta_dist_maxwell(1,idx_maxwell));

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,idx_maxwell),iii,beta_dist_maxwell(iii*2+1,idx_maxwell));
                    end
                    if fluidSetting == 'y'
                        fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                else
                    fprintf(fid,'No elastic response used.\r\n');

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,idx_maxwell),iii,beta_dist_maxwell(iii*2,idx_maxwell));
                    end
                    if fluidSetting == 'y'
                        fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                end

                % Print alternate parameters
                fprintf(fid,'\r\n\r\nAlternative Fitting Parameters (2:10)\r\n\r\n');
                for alti = 2:n_ranks
                    fprintf(fid,' \r\n\r\nRank Number %d:',alti);

                    % Print GKV
                    fprintf(fid,'\r\n\r\nGeneralized Kelvin-Voigt\r\n');
                    if elasticSetting == 'y'
                        fprintf(fid,'J_g = %3.6g [1/Pa]\r\n',beta_dist(1,I(alti)));

                        for iii = 1:n_terms
                            fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2,I(alti)),iii,beta_dist(iii*2+1,I(alti)));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,I(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist(iii*2-1,I(alti)),iii,beta_dist(iii*2,I(alti)));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist(end,I(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                    % Print Maxwell
                    fprintf(fid,'\r\n\r\nGeneralized Maxwell\r\n');
                    if elasticSetting == 'y'
                        fprintf(fid,'E_e =%3.6g [Pa]\r\n',beta_dist_maxwell(1,I_maxwell(alti)));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2+1,I_maxwell(alti)));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2,I_maxwell(alti)));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'eta_f =%3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                end

            otherwise
                error('The fit method provided is not implemented. Please check the value of fitMethod.')

        end % End switch fitMethod
        
        % Give the user some quick visual feedback showing the final fit
        if exist('dataPlot','var')
            figure(dataPlot);
            clf;
        else
            dataPlot = figure('position',[25 25 1000 1000]);
        end
        hold on
        grid on
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        set(gca,'Box','on')
        set(findall(gcf,'-property','FontSize'),'FontSize',18)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',3)
        xlabel('Time [s]', 'FontSize', 32, 'Interpreter', 'Latex')
        ylabel('Normalized Data [Various]', 'FontSize', 26)
        set(gca, 'LineStyleOrder', {'-', ':', '--', '-.'}); % different line styles
        switch avAnalysis
            case 'i'
                title(sprintf('Viscoelastic Model Performance, File No. %d',i), 'FontSize', 26)
            case 'a'
                title(sprintf('Viscoelastic Model Performance, v_{app} = %1.3g [nm/s]',v_unique(i)), 'FontSize', 26)
        end
        switch(lower(fitMethod))
            case 'iterative'
                scatter(x_fit_loop(:,1),h_norm_loop,150,'rx','linewidth',4)
                plot(x_fit_loop(:,1),F_conv_loop(dataStruct(indShift+i_loop).convParams,x_fit_loop),'k-','linewidth',3)
                scatter(maxwellInputs_loop(:,1),F_norm_loop,150,'bx','linewidth',4)
                plot(maxwellInputs_loop(:,1),h_conv_maxwell_loop(dataStruct(indShift+i_loop).maxwellParams,maxwellInputs_loop),'k--','linewidth',3)
            case 'open'
                scatter(x_fit(:,1),dataStruct(indShift+i_loop).h_norm,150,'rx','linewidth',4)
                plot(x_fit(:,1),F_conv(dataStruct(indShift+i_loop).convParams,x_fit),'k-','linewidth',3)
                scatter(maxwellInputs(:,1),dataStruct(indShift+i_loop).F_norm,150,'bx','linewidth',4)
                plot(maxwellInputs(:,1),h_conv_maxwell(dataStruct(indShift+i_loop).maxwellParams,maxwellInputs),'k--','linewidth',3)
        end
        legend({'Voigt Normalized Data',sprintf('Voigt Model, %d Term(s)',n_terms),'Maxwell Normalized Data',sprintf('Maxwell Model, %d Term(s)',n_terms)},'orientation','vertical','location','best');
        hold off

        if timeOperation
            fitSwitchEndTime(iSims,i_loop) = toc;
        end

        if strcmp(avAnalysis,'i')
            save([pathname sprintf('/%s-dataStruct-%d_terms-%s.mat',savePrepend,n_terms,saveLabel)],'dataStruct')
            fprintf('\nSaved the dataStruct to an m-file for %d terms on File No. %d\n',n_terms,i_loop);
        end

        fprintf(fid,'\r\n\r\n=================\r\n\r\n');

    end % End loopInd
            
    if timeOperation
        loopEndTimeStore(iSims) = toc;
    end
    
    fclose('all');
    
    if strcmp(avAnalysis,'a')
        save([pathname sprintf('/%s-dataStruct-%d_terms-average.mat',savePrepend,n_terms)],'dataStruct')
        fprintf('\nSaved the dataStruct to an m-file for %d terms and all velocities.\n',n_terms);
    end
            
end % End iVoigt loop

%% Save timing variables to their own file
if timeOperation
    save([pathname sprintf('/%s-TimerData-%d_terms-%s.mat',savePrepend,n_terms,avAnalysis)],...
        'fitSwitchStartTime', 'fitSwitchEndTime', 'elasticParallelStartTime', 'elasticParallelEndTime',...
        'elasticParallelBytes', 'iterativeFitStartTime', 'iterativeFitEndTime',...
        'parallelFitStartTime', 'parallelFitEndTime', 'parallelFitParallelBytes',...
        'loopStartTimeStore', 'loopEndTimeStore', 'parallelFitBreakdown',...
        'parallelFitMemoryMonitor')
    fprintf('\nSaved the timer data to an m-file.\n');
end

%% Delete Parallel Pool of MATLAB Workers
poolobj = gcp('nocreate');
delete(poolobj);