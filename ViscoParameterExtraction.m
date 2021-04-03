format long
clc

% ======================== %
% Viscoelastic Data Analysis
% Created by cparvini
% ======================== %
% This file will perform the Visco
% -elastic data analysis on 
% AFM SFS curves to extract 
% (hopefully) meaningful 
% information about how stiff
% the samples are on specific
% timescales.
% ======================== %

%% Load IGOR File and Parse
addpath Igor2Matlab
addpath(genpath([pwd '/Bruker2Matlab']))
addpath(genpath([pwd '/parTicToc']))
addpath Data

if ~exist('N_workers','var')
    N_workers = 8;
end

% Manual Programmatic settings
removeNegatives = 1;
bDisplay = 0;

% set(0,'DefaultFigureWindowStyle','normal') % Plot figures normally
% set(0,'DefaultFigureWindowStyle','docked') % Plot figures docked

filterType = 'FIR';
findRep = 'forward';
N = 2;                                  % Second order for Butterworth
cutoff_Hz = 5e3;                        % Butterworth cutoff

% Manual Simulation Setting
r_tip = 12.5e-6; % m
useBruker = 0;
n_terms_tot = 5;
n_samples = 1000;
n_maxIterations = 1e4;
multiSim = 'y';
skipNums = 0;
elasticSetting = 'y';
elasticSetting_maxwell = 'y'; % Add control to GUI (To-Do)
fluidSetting = 'n';
fluidSetting_maxwell = 'n'; % This arises from the fact that rheodictic behavior for the maxwell doesn't involve an extra fluidity term
avAnalysis = 'a';
smoothData = 'y';
scaling = 1e-60;
ignoreWarns = 1;
advancedLog = 0;
nu_sample = 0.5;
plotFigs = 0;
forwardFitTimescale = 1;
minTimescale = 1e-4;
relaxationFactor = 10;
ub_relaxationFactor = relaxationFactor;
lb_relaxationFactor = 1/(relaxationFactor);
timeOperation = 1;
fitMethod = 'iterative';
fitLog = 0;
multiLoadFitting = 0;
tossThresh = 0;
clipData = 1;
includeRetract = 0;
guiXscale = 'log';
guiYscale = 'log';
if (includeRetract) clipData = 0; end % You can't use the retract if you clip!!
weightedFit = 1;
scaleParams = 0;
paramScale = 'logscale'; % Valid Options: 'equal', 'individual', 'logscale'
optimizationMethod = 'nls'; % Valid Options: 'gd','gen-alg','nls','nelder-mead','particle-swarm','pattern','surrogate'
if ~fitLog && ~strcmp(optimizationMethod,'nls') fitLog = 1; end % The other methods already use non-log data. Leave the setting as-is
symbolicJacobian = 1;
fitWithLSQ = 1;
fitSecondDatastream = 0;
if fitSecondDatastream && fitLog
    error("You can't fit using logarithmic scale data AND use the second datastream. Please select one or the other.");
end
bruteForceBeta = 0;
forceResidual = 0;
useJacobian = 0;
fitElasticFirst = 1;
enforceGridGuess = 0;
perturbGD = 0;
normalizeGD = 0; % For using a factor of 1/N in front of the residual/sse calculation in the GD
settleGD = ''; % Choose the way you would like to settle the GD params: '', 'newton', 'nelder-mead', 'gen-alg'
saveGDdata = 1;
multiGridGD = 0;
multiGrid_hardCap = 5;
precisionGD = 1e-6;            % The threshold below which our models are "accurate"
precisionMultigrid = 1e-6;            % The threshold below which our models are "accurate"

% GD Settings
GDmethod = 'alrm';           % Options: 'standard', 'sgd', 'alrm', 'alrm-norm', 'adadelta'
alpha = 1;                  % Initial Learning Speed
alphaV = 1;                 % Initial Learning Speed
alphaS = 1;                 % Initial Learning Speed
gamma = 0.30;               % Momentum Coefficient
gammaV = 0.30;              % Momentum Coefficient
gammaS = 0.30;              % Momentum Coefficient
zeta = 1.04;                % Fraction increase beyond which the learning rate is backed off
rho = 0.75;                 % Fraction of the learning rate to use for the next step after backing off
rho_ada = 0.95;             % Decay Constant for Adadelta
rho_adaV = 0.95;            % Decay Constant for Adadelta
rho_adaS = 0.95;            % Decay Constant for Adadelta
epsilon_adadelta = 1e-6;    % Denominator conditioning parameter for Adadelta
eta = 1.15;                 % Fraction increase for the Learning Parameter if the error drops
accuracy = 1e-5;            % How closely the parameters can change
abandonThreshold = 0.90;    % If the error hasn't decreased at least X% by 1000 iterations, abandon this beta0
perturbMethod = 'LR';       % Choose quantity to perturb, 'params', or 'LR'
perturbWindow = 1000;
perturbBand = 0.05;

% Save the settings for later
settingsStruct = struct;
settingsStruct.removeNegatives = removeNegatives;
settingsStruct.bDisplay = bDisplay;
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
settingsStruct.ignoreWarns = ignoreWarns;
settingsStruct.advancedLog = advancedLog;
settingsStruct.nu_sample = nu_sample;
settingsStruct.forwardFitTimescale = forwardFitTimescale;
settingsStruct.minTimescale = minTimescale;
settingsStruct.plotFigs = plotFigs;
settingsStruct.relaxationFactor = relaxationFactor;
settingsStruct.timeOperation = timeOperation;
settingsStruct.fitMethod = fitMethod;
settingsStruct.multiLoadFitting = multiLoadFitting;
settingsStruct.optimizationMethod = optimizationMethod;
settingsStruct.fitWithLSQ = fitWithLSQ;
settingsStruct.clipData = clipData;
settingsStruct.weightedFit = weightedFit;
settingsStruct.fitSecondDatastream = fitSecondDatastream;
settingsStruct.bruteForceBeta = bruteForceBeta;
settingsStruct.forceResidual = forceResidual;
settingsStruct.useJacobian = useJacobian;
settingsStruct.fitElasticFirst = fitElasticFirst;
settingsStruct.enforceGridGuess = enforceGridGuess;
settingsStruct.GDmethod = GDmethod;
settingsStruct.scaleParams = scaleParams;
settingsStruct.paramScale_modulus = paramScale;
settingsStruct.perturbGD = perturbGD;
settingsStruct.normalizeGD = normalizeGD;
settingsStruct.settleGD = settleGD;
settingsStruct.saveGDdata = saveGDdata;
settingsStruct.multiGridGD = multiGridGD;
settingsStruct.alpha = alpha;
settingsStruct.alphaV = alphaV;
settingsStruct.alphaS = alphaS;
settingsStruct.gamma = gamma;
settingsStruct.gammaV = gammaV;
settingsStruct.gammaS = gammaS;
settingsStruct.zeta = zeta;
settingsStruct.rho = rho;
settingsStruct.rho_ada = rho_ada;
settingsStruct.rho_adaV = rho_adaV;
settingsStruct.rho_adaS = rho_adaS;
settingsStruct.eta = eta;
settingsStruct.accuracy = accuracy;
settingsStruct.precisionGD = precisionGD;
settingsStruct.abandonThreshold = abandonThreshold;
save([pathname sprintf('/%s-settingsStruct-average.mat',savePrepend)],'settingsStruct')
fprintf('\nSaved the simulation settings to an m-file.\n');
clearvars settingsStruct

FilesCheck=dir([pathname '/*.*']);

% Remove Directories
FilesCheck=FilesCheck(~ismember({FilesCheck.name},{'.','..'}));
toRemove = find([FilesCheck.isdir] == 1);
FilesCheck(toRemove) = [];

% Remove Filetypes To Ignore
toRemove = find(~endsWith({FilesCheck.name}, {'.ibw','.txt','.spm','.mat','.csv'}));
FilesCheck(toRemove) = [];

% Remove Previous Results from the count
toRemove = find(contains({FilesCheck.name}, {'FitParams'}));
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
            
        case lower('ELopez')
            Files=dir([pathname '/*.ibw']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            end
            
        case lower('ELopezSim')
            Files=dir([pathname '/*.txt']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            end
            
        case lower({'CytoDcell', '1205LuWT'})
            if strcmp(FilesTemp{end},'txt')
                Files=dir([pathname '/*Time.txt']);
            elseif strcmp(FilesTemp{end},'spm')
                Files=dir([pathname '/*.spm']);
            end
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
                
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
            
        case lower({'Control', 'treated'})
            Files=dir([pathname '/*.ibw']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{2},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{4},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            end
            
        case lower({'IterativePaper'})
            Files=dir([pathname '/*.mat']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings'}));
            Files(toRemove) = [];
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

                v_approach = str2num(FileInfo{5})*1E-9;                     % Approach Velocity, nm/s
                point_number(k) = k;                                        % Position Number
                run_number(k) = 1;                                          % Run Number
            end
         
        case lower({'KeidarCellsControl', 'KeidarCellsTreated'})
            Files=dir([pathname '/*.ibw']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{5},'-','.'))*1E-9;  % Approach Velocity, nm/s
                cell_number(k) = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
                point_number(k) = str2num(strrep(FileInfo{3},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{4},'-','.'))+1;     % Run Number
                
            end
            
        case lower({'ExcelDebugVoigt','ExcelDebugMaxwell'})
            Files=dir([pathname '/*.csv']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            end
            
        case lower({'TestCondition','TipRadiusTestCondition'})
            Files=dir([pathname '/*.mat']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings'}));
            Files(toRemove) = [];
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(FileInfo{5})*1E-9;                     % Approach Velocity, nm/s
                point_number(k) = k;                                        % Position Number
                run_number(k) = 1;                                          % Run Number
            end
            
        case lower({'NSF'})
            Files=dir([pathname '/*.csv']);
            for k=1:length(Files)
                FileNames = Files(k).name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

                v_approach(k) = str2num(strrep(FileInfo{3},'-','.'))*1E-9;  % Approach Velocity, nm/s
                point_number(k) = str2num(strrep(FileInfo{4},'-','.'));     % Position Number
                run_number(k) = str2num(strrep(FileInfo{5},'-','.'))+1;     % Run Number
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

        case lower('ELopez')
            Files=dir([pathname '/*.ibw']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number

        case lower('ELopezSim')
            Files=dir([pathname '/*.txt']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            
        case lower({'CytoDcell', '1205LuWT'})
            if strcmp(FilesTemp{end},'txt')
                Files=dir([pathname '/*Time.txt']);
                FileNames = Files.name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
            elseif strcmp(FilesTemp{end},'spm')
                Files=dir([pathname '/*.spm']);
                FileNames = Files.name;
                FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);
            end

            point_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'));     % Run Number
        
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
            
        case lower({'Control', 'treated'})
            Files=dir([pathname '/*.ibw']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{2},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{4},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
        
        case lower({'IterativePaper'})
            Files=dir([pathname '/*.mat']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings'}));
            Files(toRemove) = [];
            
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

            v_approach = str2num(FileInfo{5})*1E-9;  % Approach Velocity, nm/s
            point_number = 1;                                        % Position Number
            run_number = 1;                                          % Run Number
           
        case lower({'KeidarCellsControl', 'KeidarCellsTreated'})
            Files=dir([pathname '/*.ibw']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{5},'-','.'))*1E-9;  % Approach Velocity, nm/s
            cell_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            point_number = str2num(strrep(FileInfo{3},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{4},'-','.'))+1;     % Run Number
            
        case lower({'ExcelDebugVoigt','ExcelDebugMaxwell'})
            Files=dir([pathname '/*.csv']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            
        case lower({'TestCondition','TipRadiusTestCondition'})
            Files=dir([pathname '/*.mat']);
            toRemove = find(contains({Files.name}, {'settingsStruct','Settings'}));
            Files(toRemove) = [];
            
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '-' '.'},'CollapseDelimiters',true);

            v_approach = str2num(FileInfo{5})*1E-9;                 % Approach Velocity, nm/s
            point_number = 1;                                       % Position Number
            run_number = 1;                                         % Run Number
            
        case lower({'NSF'})
            Files=dir([pathname '/*.csv']);
            FileNames = Files.name;
            FileInfo = strsplit(FileNames, {'_' '.'},'CollapseDelimiters',true);

            v_approach = str2num(strrep(FileInfo{3},'-','.'))*1E-9;  % Approach Velocity, nm/s
            point_number = str2num(strrep(FileInfo{4},'-','.'));     % Position Number
            run_number = str2num(strrep(FileInfo{5},'-','.'))+1;     % Run Number
            
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
%             k_cantilever(k) = 40.73; % Value in N/m
            
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
    
        case lower('ELopez')
            RawData = IBWread([Files(k).folder '/' Files(k).name]);
            [headerValue,null] = strsplit(RawData.WaveNotes,'\r',...
            'DelimiterType','RegularExpression');

            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, m/s
            point_number_temp = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number_temp = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            
            dataStruct(k).z = RawData.y(:,5);
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
                        
        case lower('ELopezSim')
            
            if length(Files) > 1
                fileID = fopen(fullfile([Files(k).folder '/' Files(k).name]), 'rt');
                RawData = textscan(fileID,'%f%f%f%f%f','Delimiter',' ','HeaderLines',1);
                fclose(fileID);
            else
                fileID = fopen([Files.folder '/' Files.name], 'rt');
                RawData = textscan(fileID,'%f%f%f%f%f','Delimiter',' ','HeaderLines',1);
                fclose(fileID);
            end
            
            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;  % Approach Velocity, m/s
            point_number_temp = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            run_number_temp = str2num(strrep(FileInfo{3},'-','.'))+1;     % Run Number
            
            dataStruct(k).z = RawData{2}.*1E-9;
            dataStruct(k).d = RawData{5}.*1E-9;
            
            k_cantilever(k) = 0.1; % Value in N/m

            % Insert Fake InVOLS (Not needed for this)
            defl_InVOLS(k) = 1; % Value in nm/V
                
            % Create time array
            dataStruct(k).dt = mean(diff(RawData{1})); % Value in s
            dataStruct(k).n_sam = size(RawData{1},2);
            dataStruct(k).t = RawData{1};
            
            manualParams = [2.0e-10, 1.0e-9, 1.0e-4, 5.0e-8, 1.0e-3, 1.0e-7, 1.0e-2, 1.0e-6, 1.0e-1, 1.0e-6, 1.0e0, 5e-16];

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

                %----------------------------
                %Force timed plot
    %             f = figure();
    %             movegui(f, 'northwest');
    % 
    %             %get timed data (Volts) of channel 1
                [time, deflectionError, xLabel, yLabel] = NSMU.CreateForceTimePlot(1, NSMU.VOLTS);
    %             [chanDesc] = NSMU.GetDataTypeDesc(1);
    %             plot(time, deflectionError);
    % 
    %             title (strcat(chanDesc, ' vs. Time'));
    %             xlabel(xLabel);
    %             ylabel(yLabel);

                %Get F vs Z plot of channel 1
                [zApproach, zRetract, FApproach, FRetract, xLabel, yLabel] = NSMU.CreateForceZPlot(1, NSMU.FORCE, 0);

                %---------------------------
                % Force z plot
    %             f = figure();
    %             movegui(f, 'north');
    %             hold on;
    %             plot(zApproach, FApproach,'r');
    %             plot(zRetract, FRetract,'b');
    %             title ('Force vs. Z');
    %             xlabel(xLabel);
    %             ylabel(yLabel);

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

                %---------------------------
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

                % nu_sample = NSMU.GetPoissonRatio();
                
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
                    
                    % See what section of data we're actually sending
                    % plot(zApproach*(1e-6))
                    % hold on
                    % plot(dataStruct(k).z,'r--')
                    % hold off
                else
                    dataStruct(k).z = zApproach(dataError_pos+buffer:end)*(1e-6);
                    dataStruct(k).d = approachDefl(dataError_pos+buffer:end)*(1e-9);
                    dataStruct(k).t = time(1:(end-dataError_pos+buffer));
                    dataStruct(k).dt = dt;
                    dataStruct(k).n_sam = length(dataStruct(k).t);

                    % See what section of data we're actually sending
                    % plot(zApproach*(1e-6))
                    % hold on
                    % plot(((dataError_pos+buffer):length(zApproach)),dataStruct(k).z,'r--')
                    % hold off
                end
                
                v_approach(k) = mean(abs(diff(dataStruct(k).z))./(dataStruct(k).dt*10));    % Approach Velocity, m/s
                
                fclose(fid);
                
            end
            
        case lower({'Control', 'treated'})
            RawData = IBWread([Files(k).folder '/' Files(k).name]);
            [headerValue,null] = strsplit(RawData.WaveNotes,'\r',...
            'DelimiterType','RegularExpression');

            v_approach_temp = str2num(strrep(FileInfo{2},'-','.'))*1E-9;  % Approach Velocity, m/s
            point_number_temp = str2num(strrep(FileInfo{4},'-','.'));     % Position Number
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
            
        case lower({'IterativePaper'})
            RawData = load([Files(k).folder '/' Files(k).name],'z','d','t','F');
            
            settingsCheck=dir([pathname '/*.*']);

            % Remove Directories
            settingsCheck=settingsCheck(~ismember({settingsCheck.name},{'.','..'}));
            toRemove = find([settingsCheck.isdir] == 1);
            settingsCheck(toRemove) = [];

            % Remove Filetypes To Ignore
            toRemove = find(~endsWith({settingsCheck.name}, {'.mat'}));
            settingsCheck(toRemove) = [];

            toRemove = find(~contains({settingsCheck.name}, {'Settings'}));
            settingsCheck(toRemove) = [];
            
            settingsData = load([settingsCheck.folder '/' settingsCheck.name]);
            
            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;    % Approach Velocity, m/s
            point_number_temp = k;                                          % Position Number
            run_number_temp = 1;                                            % Run Number

            dataStruct(k).z = -(RawData.z-RawData.z(1));
            dataStruct(k).t = RawData.t;
            dataStruct(k).n_sam = numel(dataStruct(k).t);
            
            k_cantilever(k) = settingsData.settingsStruct.k_m1;
            dataStruct(k).d = RawData.F./k_cantilever(k);
            defl_InVOLS(k) = 1;
            dataStruct(k).dt = settingsData.settingsStruct.dt;
            
        case lower({'KeidarCellsControl', 'KeidarCellsTreated'})
            RawData = IBWread([Files(k).folder '/' Files(k).name]);
            [headerValue,null] = strsplit(RawData.WaveNotes,'\r',...
            'DelimiterType','RegularExpression');

            v_approach_temp = str2num(strrep(FileInfo{5},'-','.'))*1E-9;  % Approach Velocity, m/s
            cell_number_temp = str2num(strrep(FileInfo{2},'-','.'));     % Position Number
            point_number_temp = str2num(strrep(FileInfo{3},'-','.'));     % Position Number
            run_number_temp = str2num(strrep(FileInfo{4},'-','.'))+1;     % Run Number
            
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

            [z_max_temp, z_max_ind_temp] = max(dataStruct(k).z);
            v_approach(k) = round(mean(diff(dataStruct(k).z(1:z_max_ind_temp))./dataStruct(k).dt),1,'significant');
            
        case lower({'ExcelDebugVoigt','ExcelDebugMaxwell'})
            opts = detectImportOptions([Files(k).folder '/' Files(k).name]);
%             opts.NumHeaderLines = 1;
%             opts.OutputType = 'double';
            RawData = readmatrix([Files(k).folder '/' Files(k).name],opts,'OutputType','double');
            RawData(isnan(RawData(:,1)),:) = [];
            dataStruct(k).dt = 2e-5;%RawData(2,1)-RawData(1,1);
%             t_temp = RawData(:,1);
%             z_temp = RawData(:,2);
%             d_temp = RawData(:,5);
            t_interp = (dataStruct(k).dt:dataStruct(k).dt:max(RawData(:,1)))';
            
            [~, uniqueZ] = unique(RawData(:,2));
            [~, uniqueD] = unique(RawData(:,5));
            
            z_interp = interp1(RawData(uniqueZ,1), RawData(uniqueZ,2),...
                    t_interp, 'linear', NaN); % Interploate To New Time Values
            d_interp = interp1(RawData(uniqueD,1), RawData(uniqueD,5),...
                    t_interp, 'linear', NaN); % Interploate To New Time Values
            
            nanCheck = isnan(reshape(t_interp,length(t_interp),1)) +...
                isnan(reshape(z_interp,length(z_interp),1)) +...
                isnan(reshape(d_interp,length(d_interp),1));
            nanCheck(nanCheck ~= 0) = 1;
            
            t_interp(nanCheck == 1) = [];
            z_interp(nanCheck == 1) = [];
            d_interp(nanCheck == 1) = [];
                
            if size(z_interp,2) > size(z_interp,1)
                z_interp = z_interp';
            end

            if size(d_interp,2) > size(d_interp,1)
                d_interp = d_interp';
            end

            if size(t_interp,2) > size(t_interp,1)
                t_interp = t_interp';
            end

            dataStruct(k).t = t_interp;
            dataStruct(k).z = -z_interp;
            dataStruct(k).d = d_interp;
            dataStruct(k).n_sam = numel(dataStruct(k).t);
            k_cantilever(k) = 2; % N/m
            defl_InVOLS(k) = 1; % Fake value
            
            manualParams = [1e6 1e5 1e-3];
            
        case lower({'TestCondition','TipRadiusTestCondition'})
            RawData = load([Files(k).folder '/' Files(k).name],'z','d','t','F');
                        
            settingsCheck=dir([pathname '/*.*']);

            % Remove Directories
            settingsCheck=settingsCheck(~ismember({settingsCheck.name},{'.','..'}));
            toRemove = find([settingsCheck.isdir] == 1);
            settingsCheck(toRemove) = [];

            % Remove Filetypes To Ignore
            toRemove = find(~endsWith({settingsCheck.name}, {'.mat'}));
            settingsCheck(toRemove) = [];

            toRemove = find(~contains({settingsCheck.name}, {'Settings'}));
            settingsCheck(toRemove) = [];
            
            if size(settingsCheck,1) > 1
                setInd = k;
            else
                setInd = 1;
            end
            
            settingsData = load([settingsCheck(setInd).folder '/' settingsCheck(setInd).name]);
            
            v_approach_temp = str2num(strrep(FileInfo{4},'-','.'))*1E-9;    % Approach Velocity, m/s
            point_number_temp = k;                                          % Position Number
            run_number_temp = 1;                                            % Run Number

            dataStruct(k).z = -(RawData.z-RawData.z(1));
            dataStruct(k).t = RawData.t;
            dataStruct(k).n_sam = numel(dataStruct(k).t);
            
            k_cantilever(k) = settingsData.settingsStruct.k_m1;
            dataStruct(k).d = RawData.F./k_cantilever(k);
            defl_InVOLS(k) = 1;
            dataStruct(k).dt = settingsData.settingsStruct.dt;
            
            condSwitch = FileInfo{1};
            condSwitch = condSwitch(isstrprop(condSwitch,'digit'));
            manualParams = {};
            switch str2num(condSwitch)
                case 1
                    manualParams{k} = [1e4 1e4 1e-4];
                case 2
                    manualParams{k} = [1e4 1e4 1e-4 1e4 1e-3];
                case 3
                    manualParams{k} = [1e4 1e4 1e-4 1e4 1e-3 1e4 1e-2];
                case 4
                    manualParams{k} = [1e4 1e4 1e-4 1e4 1e-3 1e4 1e-2 1e4 1e-1];
                case 5
                    manualParams{k} = [5e5 1e4 1e-4];
                case 6
                    manualParams{k} = [5e5 1e4 1e-4 1e4 1e-3];
                case 7
                    manualParams{k} = [5e5 1e4 1e-4 1e4 1e-3 1e4 1e-2];
                case 8
                    manualParams{k} = [5e5 1e4 1e-4 1e4 1e-3 1e4 1e-2 1e4 1e-1];
                case 9
                    manualParams{k} = [5e7 1e4 1e-4];
                case 10
                    manualParams{k} = [5e7 1e4 1e-4 1e4 1e-3];
                case 11
                    manualParams{k} = [5e7 1e4 1e-4 1e4 1e-3 1e4 1e-2];
                case 12
                    manualParams{k} = [5e7 1e4 1e-4 1e4 1e-3 1e4 1e-2 1e4 1e-1];
                case 13
                    manualParams{k} = [1e4 5e5 1e-4];
                case 14
                    manualParams{k} = [1e4 5e5 1e-4 5e5 1e-3];
                case 15
                    manualParams{k} = [1e4 5e5 1e-4 5e5 1e-3 5e5 1e-2];
                case 16
                    manualParams{k} = [1e4 5e5 1e-4 5e5 1e-3 5e5 1e-2 5e5 1e-1];
                case 17
                    manualParams{k} = [5e5 5e5 1e-4];
                case 18
                    manualParams{k} = [5e5 5e5 1e-4 5e5 1e-3];
                case 19
                    manualParams{k} = [5e5 5e5 1e-4 5e5 1e-3 5e5 1e-2];
                case 20
                    manualParams{k} = [5e5 5e5 1e-4 5e5 1e-3 5e5 1e-2 5e5 1e-1];
                case 21
                    manualParams{k} = [5e7 5e5 1e-4];
                case 22
                    manualParams{k} = [5e7 5e5 1e-4 5e5 1e-3];
                case 23
                    manualParams{k} = [5e7 5e5 1e-4 5e5 1e-3 5e5 1e-2];
                case 24
                    manualParams{k} = [5e7 5e5 1e-4 5e5 1e-3 5e5 1e-2 5e5 1e-1];
                case 25
                    manualParams{k} = [1e4 5e7 1e-4];
                case 26
                    manualParams{k} = [1e4 5e7 1e-4 5e7 1e-3];
                case 27
                    manualParams{k} = [1e4 5e7 1e-4 5e7 1e-3 5e7 1e-2];
                case 28
                    manualParams{k} = [1e4 5e7 1e-4 5e7 1e-3 5e7 1e-2 5e7 1e-1];
                case 29
                    manualParams{k} = [5e5 5e7 1e-4];
                case 30
                    manualParams{k} = [5e5 5e7 1e-4 5e7 1e-3];
                case 31
                    manualParams{k} = [5e5 5e7 1e-4 5e7 1e-3 5e7 1e-2];
                case 32
                    manualParams{k} = [5e5 5e7 1e-4 5e7 1e-3 5e7 1e-2 5e7 1e-1];
                case 33
                    manualParams{k} = [5e7 5e7 1e-4];
                case 34
                    manualParams{k} = [5e7 5e7 1e-4 5e7 1e-3];
                case 35
                    manualParams{k} = [5e7 5e7 1e-4 5e7 1e-3 5e7 1e-2];
                case 36
                    manualParams{k} = [5e7 5e7 1e-4 5e7 1e-3 5e7 1e-2 5e7 1e-1];
            end
            
        case lower({'NSF'})
            opts = detectImportOptions([Files(k).folder '/' Files(k).name]);
%             opts.NumHeaderLines = 1;
%             opts.OutputType = 'double';
            RawData = readmatrix([Files(k).folder '/' Files(k).name],opts,'OutputType','double');
%             RawData(isnan(RawData(:,1)),:) = [];
            dataStruct(k).dt = RawData(2,1)-RawData(1,1);

            dataStruct(k).t = RawData(:,1);
            dataStruct(k).z = RawData(:,2);
            [~, z_max_ind] = max(-dataStruct(k).z);
            dataStruct(k).z_approach = dataStruct(k).z(1:z_max_ind);
            
%             repInd = find((dataStruct(k).z)<0,1);
            repInd = (find(flipud(RawData(:,3))<0,1) - 1);
            if ~isempty(repInd)
                repInd = length(RawData(:,3)) - repInd + 1;
                dataStruct(k).h_r = -(RawData(repInd:end,2)-RawData(repInd-1,2));
            else
                repInd = 1;
                dataStruct(k).h_r = -(RawData(repInd:end,2)-RawData(repInd,2));
            end
            
            dataStruct(k).F_r = RawData(repInd:end,3);
            dataStruct(k).F_r_smooth = dataStruct(k).F_r;
            dataStruct(k).t_r = ( (1:length(RawData(repInd:end,1))).*dataStruct(k).dt )';
            dataStruct(k).h_r_smooth = dataStruct(k).h_r;
            dataStruct(k).n_sam = numel(dataStruct(k).t);
            
            % Change to Logarithmic Scaling to use Enrique et al.'s Fit Method
            tr = dataStruct(k).dt;
            st = dataStruct(k).t(end);

            F_log = [];
            t_log = [];

            if ~advancedLog
                F_log = log_scale(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
                t_log = log_scale(dataStruct(k).t_r,dataStruct(k).t_r,tr,st);
            else
                [F_log, ind_set]= log_scale_advanced(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
                t_log = dataStruct(k).t_r(ind_set)';
            end

            % Save the Log-Form of the tip force (F) and time array (t_r)
            dataStruct(k).F_r_log = F_log;
            dataStruct(k).t_r_log = t_log;
                        
            dataStruct(k).z_corrected = RawData(:,2);
            dataStruct(k).d_corrected = zeros(size(RawData(:,2)));
            dataStruct(k).t_approach = RawData(:,1);
            
            % Copy to averaged row
            dataStruct(k+1).dt = dataStruct(k).dt;
            dataStruct(k+1).t = dataStruct(k).t;
            dataStruct(k+1).z = dataStruct(k).z;
            dataStruct(k+1).z_approach = dataStruct(k).z_approach;
            dataStruct(k+1).t_r = dataStruct(k).t_r;
            dataStruct(k+1).h_r = dataStruct(k).h_r;
            dataStruct(k+1).h_r_smooth = dataStruct(k).h_r_smooth;
            dataStruct(k+1).F_r = dataStruct(k).F_r;
            dataStruct(k+1).F_r_smooth = dataStruct(k).F_r_smooth;
            dataStruct(k+1).n_sam = dataStruct(k).n_sam;
            dataStruct(k+1).F_r_log = dataStruct(k).F_r_log;
            dataStruct(k+1).t_r_log = dataStruct(k).t_r_log;
            dataStruct(k+1).z_corrected = dataStruct(k).z_corrected;
            dataStruct(k+1).d_corrected = dataStruct(k).d_corrected;
            dataStruct(k+1).t_approach = dataStruct(k).t_approach;
                        
            continue;
            
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
        
        % Check out the filter performance (frequency domain)
%         f = linspace(0,1/(dataStruct(k).dt),length(dataStruct(k).d_corrected));
%         figure
%         loglog(f,abs(fft(dataStruct(k).d_corrected)))
%         hold on
%         loglog(f,abs(fft(LPF(dataStruct(k).d_corrected))))
%         grid on
%         hold off
        
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
    
    if ~advancedLog
        F_log = log_scale(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
        t_log = log_scale(dataStruct(k).t_r,dataStruct(k).t_r,tr,st);
    else
        [F_log, ind_set]= log_scale_advanced(dataStruct(k).F_r,dataStruct(k).t_r,tr,st);
        t_log = dataStruct(k).t_r(ind_set)';
    end
    
    % Save the Log-Form of the tip force (F) and time array (t_r)
    dataStruct(k).F_r_log = F_log;
    dataStruct(k).t_r_log = t_log;
    
%     clearvars -except k dataStruct tip_rep_pos_all tip_rep_pos_all_smooth...
%         dSmoothMinAll Files searchpath pathname savePrepend removeNegatives...
%         bDisplay filterType findRep N cutoff_Hz r_tip element_type...
%         n_terms_tot n_terms_maxwell n_samples n_maxIterations fitDataSource...
%         multiSim skipNums elasticSetting fluidSetting avAnalysis smoothData...
%         scaling saveMats ignoreWarns advancedLog run_number point_number...
%         v_approach k_cantilever point_number_temp run_number_temp useBruker...
%         nu_sample app plotFigs relaxationFactor fitHarmonics forwardFitTimescale...
%         ub_relaxationFactor lb_relaxationFactor timeOperation

end

clearvars k point_number_temp run_number_temp

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
    
    if strcmp(strSwitch,'NSF')
       break;
    end
    
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
    
%     dataStruct(length(Files)+i).z_approach = dataStruct(length(Files)+i).z_average(1:z_max_ind);
%     dataStruct(length(Files)+i).d_approach = dataStruct(length(Files)+i).d_average(1:z_max_ind);
    
    % Save the data for later, if you want
%     dt = dataStruct(length(Files)+i).dt;
%     t = dataStruct(length(Files)+i).t_average;
%     z = dataStruct(length(Files)+i).z_average;
%     d = dataStruct(length(Files)+i).d_average;
%     save([pathname sprintf('/%s-Data-%.f_nm-s.mat',savePrepend,v_unique(i)/1E-9)],'dt','t','z','d');
    
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
        
        % Check out the filter performance (frequency domain)
%         f = linspace(0,1/(dataStruct(length(Files)+i).dt),length(dataStruct(length(Files)+i).d_approach));
%         figure
%         loglog(f,abs(fft(dataStruct(length(Files)+i).d_approach)))
%         hold on
%         loglog(f,abs(fft(LPF(dataStruct(length(Files)+i).d_approach))))
%         grid on
%         hold off
        
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

    if ~advancedLog
        F_log = log_scale(dataStruct(length(Files)+i).F_r,dataStruct(length(Files)+i).t_r,tr,st);
        t_log = log_scale(dataStruct(length(Files)+i).t_r,dataStruct(length(Files)+i).t_r,tr,st);
    else
        [F_log, ind_set]= log_scale_advanced(dataStruct(length(Files)+i).F_r,dataStruct(length(Files)+i).t_r,tr,st);
        t_log = dataStruct(length(Files)+i).t_r(ind_set)';
    end
    
    % Save the Log-Form of the tip force (F) and time array (t_r)
    dataStruct(length(Files)+i).F_r_log = F_log;
    dataStruct(length(Files)+i).t_r_log = t_log;

end

clearvars d_r_average t_r_average z_r_average N z_max_ind_temp temp d_temp t_temp z_temp xi yi

%% Visualize the Data Together
close all

% legendString = {};
% figure(1)
% clf
% set(gcf, 'Position', get(0, 'Screensize'));
% hold on
% grid on
% 
% titleString = sprintf('Treated Deflection vs. Z-Sensor\nv_{app} = %4.2f [nm/s]',1E9*v_unique(i));
% 
% if length(Files) > 1
%     for ii = 1:length(Files)
%         if v_approach(ii) == v_unique(i)
%             legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number(ii),run_number(ii))});
%             plot(dataStruct(ii).z_full_corrected,dataStruct(ii).d_full_corrected)
%         else
%             continue;
%         end
%     end
% else
%     legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number,run_number)});
%     plot(dataStruct(i).z_full_corrected,dataStruct(i).d_full_corrected)
% end
% 
% % legendString = vertcat(legendString, {sprintf('Averaged Curve')});
% % plot(dataStruct(length(Files)+i).z_full_average,dataStruct(length(Files)+i).d_full_average,'k--','LineWidth',5)
% 
% set(gca, 'xlim', [-1e-7 0.55e-7])
% set(findall(gcf,'-property','FontSize'),'FontSize',28)
% set(gca,'TickLength',[0.02 0.02],'LineWidth',3)
% xlabel('Z-sensor [m]', 'FontSize', 32)
% ylabel('Deflection [m]', 'FontSize', 32)
% title(titleString, 'FontSize', 32)
% % legend(legendString,'Location','eastoutside','Orientation','Vertical', 'FontSize', 10)
% hold off
% 
% % saveas(gcf, [pathname sprintf('/AFM-SFS-CorrectedDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
% % savefig(gcf, [pathname sprintf('/AFM-SFS-CorrectedDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')

if plotFigs
    for i = 1:length(v_unique)

        legendString = {};
        figure(1)
        clf
        set(gcf, 'Position', get(0, 'Screensize'));
        hold on
        grid on

        titleString = sprintf('Deflection vs. Z-Sensor\nv_{app} = %4.2f [nm/s]',1E9*v_unique(i));

        if length(Files) > 1
            for ii = 1:length(Files)
                if v_approach(ii) == v_unique(i)
                    legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number(ii),run_number(ii))});
                    plot(dataStruct(ii).z_corrected,dataStruct(ii).d_corrected)
                else
                    continue;
                end
            end
        else
            legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number,run_number)});
            plot(dataStruct(i).z_corrected,dataStruct(i).d_corrected)
        end

        legendString = vertcat(legendString, {sprintf('Averaged Curve')});
        plot(dataStruct(length(Files)+i).z_average,dataStruct(length(Files)+i).d_average,'k--','LineWidth',5)

        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',3)
        xlabel('Z-sensor [m]', 'FontSize', 32)
        ylabel('Deflection [m]', 'FontSize', 32)
        title(titleString, 'FontSize', 32)
%         legend(legendString,'Location','eastoutside','Orientation','Vertical', 'FontSize', 10)
        hold off

%         saveas(gcf, [pathname sprintf('/AFM-SFS-CorrectedDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
%         savefig(gcf, [pathname sprintf('/AFM-SFS-CorrectedDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end

    clearvars legendString titleString
    % close all

    for i = 1:length(v_unique)

        legendString = {};
        figure(2)
        clf
    %     figure('units','normalized','outerposition',[0 0 1 1])
    %     set(gcf, 'Position',  [10, 50, 1000, 1000])
        hold on
        grid on

        titleString = sprintf('Repulsive Deflection vs. Repulsive Z-Sensor\nv_{app} = %4.2f [nm/s]',1E9*v_unique(i));

        if length(Files) > 1
            for ii = 1:length(Files)
                if v_approach(ii) == v_unique(i)
                    legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number(ii),run_number(ii))});
                    if strcmp(smoothData,'y')
                        plot(dataStruct(ii).z_r_smooth,dataStruct(ii).d_r_smooth,'LineWidth',5)
                    else
                        plot(dataStruct(ii).z_r,dataStruct(ii).d_r,'LineWidth',5)
                    end
    %                 fprintf('Plotted File %d\n',ii)
                else
    %                 fprintf('Skipped File %d\n',ii)
                    continue;
                end
            end
        else 
            legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number,run_number)});
            if strcmp(smoothData,'y')
                plot(dataStruct(i).z_r_smooth,dataStruct(i).d_r_smooth,'LineWidth',5)
            else
                plot(dataStruct(i).z_r,dataStruct(i).d_r,'LineWidth',5)
            end
        end

        legendString = vertcat(legendString, {sprintf('Averaged Curve')});
        if strcmp(smoothData,'y')
            plot(dataStruct(length(Files)+i).z_r_smooth,dataStruct(length(Files)+i).d_r_smooth,'k--','LineWidth',5)
        else
            plot(dataStruct(length(Files)+i).z_r,dataStruct(length(Files)+i).d_r,'k--','LineWidth',5)
        end

        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',3)
        xlabel('Repulsive Z-sensor [m]', 'FontSize', 32)
        ylabel('Repulsive Deflection [m]', 'FontSize', 32)
        title(titleString, 'FontSize', 20)
        legend(legendString,'Location','eastoutside','Orientation','Vertical', 'FontSize', 10)
        hold off

        saveas(gcf, [pathname sprintf('/AFM-SFS-RepulsiveDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [pathname sprintf('/AFM-SFS-RepulsiveDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end

    for i = 1:length(v_unique)

        legendString = {};
        figure(3)
        clf
    %     figure('units','normalized','outerposition',[0 0 1 1])
    %     set(gcf, 'Position',  [10, 50, 1000, 1000])
        hold on
        grid on

        titleString = sprintf('Repulsive Indentation vs. Time\nv_{app} = %4.2f [nm/s]',1E9*v_unique(i));

        if length(Files) > 1
            for ii = 1:length(Files)
                if v_approach(ii) == v_unique(i)
                    legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number(ii),run_number(ii))});
                    if strcmp(smoothData,'y')
                        plot(dataStruct(ii).t_r_smooth,dataStruct(ii).h_r_smooth, 'LineWidth', 5);
                    else
                        plot(dataStruct(ii).t_r,dataStruct(ii).h_r, 'LineWidth', 5);
                    end
    %                 fprintf('Plotted File %d\n',ii)
                else
    %                 fprintf('Skipped File %d\n',ii)
                    continue;
                end
            end
        else 
            legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number,run_number)});
            if strcmp(smoothData,'y')
                plot(dataStruct(i).t_r_smooth,dataStruct(i).h_r_smooth)
            else
                plot(dataStruct(i).t_r,dataStruct(i).h_r)
            end
        end

        legendString = vertcat(legendString, {sprintf('Averaged Curve')});
        if strcmp(smoothData,'y')
            plot(dataStruct(length(Files)+i).t_r_smooth,dataStruct(length(Files)+i).h_r_smooth,'k--','LineWidth',5)
        else
            plot(dataStruct(length(Files)+i).t_r,dataStruct(length(Files)+i).h_r,'k--','LineWidth',5)
        end

        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',3)
        xlabel('Repulsive Time [s]', 'FontSize', 32)
        ylabel('Repulsive Indentation [m]', 'FontSize', 32)
        title(titleString, 'FontSize', 32)
        legend(legendString,'Location','eastoutside','Orientation','Vertical', 'FontSize', 10)
        hold off

        saveas(gcf, [pathname sprintf('/AFM-SFS-IndentationDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [pathname sprintf('/AFM-SFS-IndentationDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end

    for i = 1:length(v_unique)

        legendString = {};
        figure(4)
        clf
    %     figure('units','normalized','outerposition',[0 0 1 1])
    %     set(gcf, 'Position',  [10, 50, 1000, 1000])
        hold on
        grid on

        titleString = sprintf('Repulsive Force vs. Time\nv_{app} = %4.2f [nm/s]',1E9*v_unique(i));

        if length(Files) > 1
            for ii = 1:length(Files)
                if v_approach(ii) == v_unique(i)
                    legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number(ii),run_number(ii))});
                    if strcmp(smoothData,'y')
                        plot(dataStruct(ii).t_r_smooth,dataStruct(ii).F_r_smooth, 'LineWidth', 5);
                    else
                        plot(dataStruct(ii).t_r,dataStruct(ii).F_r, 'LineWidth', 5);
                    end
    %                 fprintf('Plotted File %d\n',ii)
                else
    %                 fprintf('Skipped File %d\n',ii)
                    continue;
                end
            end
        else 
            legendString = vertcat(legendString, {sprintf('Point %d, Run %d',point_number,run_number)});
            if strcmp(smoothData,'y')
                plot(dataStruct(i).t_r_smooth,dataStruct(i).F_r_smooth, 'LineWidth', 5);
            else
                plot(dataStruct(i).t_r,dataStruct(i).F_r, 'LineWidth', 5);
            end
        end

        legendString = vertcat(legendString, {sprintf('Averaged Curve')});
        if strcmp(smoothData,'y')
            plot(dataStruct(length(Files)+i).t_r_smooth,dataStruct(length(Files)+i).F_r_smooth,'k--','LineWidth',5)
        else
            plot(dataStruct(length(Files)+i).t_r,dataStruct(length(Files)+i).F_r,'k--','LineWidth',5)
        end

        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        set(gca,'TickLength',[0.02 0.02],'LineWidth',3)
        xlabel('Repulsive Time [s]', 'FontSize', 32)
        ylabel('Repulsive Force [N]', 'FontSize', 32)
        title(titleString, 'FontSize', 32)
        legend(legendString,'Location','eastoutside','Orientation','Vertical', 'FontSize', 10)
        hold off

        saveas(gcf, [pathname sprintf('/AFM-SFS-ForceDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [pathname sprintf('/AFM-SFS-ForceDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end
end

%% Initialize
% for iii = 1:length(Files)
%    if dataStruct(iii).h_r(end) < 0
%       disp(iii); 
%    end
% end

if multiSim == 'y'
    nSims = n_terms_tot;
    voigtArray = 1:n_terms_tot;
else
    nSims = 1;
    voigtArray = n_terms_tot;
    skipNums = 0;
end

% Area Normalization Settings
% Using data for Silicon:
nu_tip = 0.27; % Incompressible = 0.5
E_tip = 160E9; % Young's Modulus of probe tip

% Check for 0 input from the GUI---means the user did not set it, and we
% did not update it during the data correction steps
if nu_sample == 0
    % For Human Skin:
    nu_sample = 0.5;
end

% Create timing Arrays
if timeOperation
    loopStartTimeStore = NaN(nSims,1);
    loopEndTimeStore = NaN(size(loopStartTimeStore));
    
    % Start the initial timer
    tic
end

%% Outer Loop
%% Outer Loop
fclose('all');

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
    parallelFitBreakdown = []; % Must be empty so it takes the correct datatype (Par object)?
%     parallelFitBreakdown = cell(nSims,numel(loopInds));
    parallelFitMemoryMonitor = cell(nSims,numel(loopInds));
end

for iVoigt = 1:nSims
    
    if timeOperation
        loopStartTimeStore(iVoigt) = toc;
    end
    
    fprintf('\nRunning Simulation for %d Elements\n\n',voigtArray(iVoigt));
    n_terms = voigtArray(iVoigt);
    
    if ismember(n_terms,skipNums)
       fprintf('\nSkipped %d Retardance Terms per user request.\n\n',voigtArray(iVoigt));
       continue; 
    end

    % Open Parallel Pool of MATLAB Workers
    if isempty(gcp('nocreate'))
        %parpool(N_workers)
        parpool('IdleTimeout', Inf)
    end

    fid = fopen([pathname sprintf('/FitParams-nBranches_%d.txt',voigtArray(iVoigt))],'w');
    fprintf(fid,'Fitting Parameters, %s\r\n=================\r\n', date);
    fprintf(fid,'\r\nTip Radius: %4.3g\r\nVoigt Terms: %d\r\nScaling Factor: %4.3g',r_tip,n_terms,scaling);

    switch avAnalysis
        case 'i'
            fprintf(fid,'\r\nIndividual File Analysis\r\n=================\r\n');
        case 'a'
            fprintf(fid,'\r\nAveraged File Analysis\r\n=================\r\n');
    end
    
    if ~multiLoadFitting
        
        for i_loop = loopInds

            % Make the save label
            switch avAnalysis
                case 'i'
                    saveLabel = sprintf('FileNo_%d',i_loop);
                case 'a'
                    saveLabel = sprintf('LoadLevel%d',i_loop);
            end

            x_lin = dataStruct(indShift+i_loop).t_r_log; % X data for linear fit
            y_lin = dataStruct(indShift+i_loop).F_r_log; % Y data for linear fit

            init_params = polyfit(x_lin, y_lin, 1); % Initial Linear Fit Parameters (guess)

            % Linear Fit Model Definition
            F_linFit = @(c,input) c .* input;

            lb = 0; % Zero Newtons/s is the minimum force rate possible
            ub = 1e4; % 1e4 Newton/s is the upper bound

            options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                'MaxFunctionEvaluations', n_maxIterations,...
                'MaxIterations', n_maxIterations,...
                'FiniteDifferenceType','central',...
                'FunctionTolerance', scaling,...
                'OptimalityTolerance', scaling,...
                'StepTolerance', scaling,...
                'Display', 'none');

            beta_dist = zeros(length(lb),n_samples);
            beta0_dist = zeros(size(beta_dist));
            resnorm_dist = zeros(1,n_samples);

            % Parallel Loop for Grid Search
            % Turned off parallel. To turn on, replace for with parfor.
            for i = 1:n_samples
                try
                    beta0 = init_params(1);
                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                        lsqcurvefit(F_linFit,beta0,x_lin,y_lin,lb,ub,options);

                    beta_dist(:,i) = betatemp;
                    beta0_dist(:,i) = beta0;
                    resnorm_dist(i) = sum(abs(residual));
                catch
                    beta_dist(:,i) = NaN;
                    beta0_dist(:,i) = NaN;
                    resnorm_dist(i) = NaN;
                end
                
            end

            [~,idx] = min(resnorm_dist);
            dataStruct(indShift+i_loop).Fdot = beta_dist(1,idx);
            
            clearvars x_lin y_lin beta_dist beta0_dist resnorm

            % Calculate the Linear Fluidity (chi)
            % According to Eq. 14, This is the relation between chi and tip
            % location when force is assumed to be linear
            dt = (dataStruct(indShift+i_loop).dt);
            st = dataStruct(indShift+i_loop).t_r(end);

            if strcmp(smoothData,'n')
                dataStruct(indShift+i_loop).chi_linear = ((16.0*sqrt(r_tip))/3)...
                    .* (dataStruct(indShift+i_loop).h_r.^1.5) ./ dataStruct(indShift+i_loop).Fdot;
                dataStruct(indShift+i_loop).chi_linear_log = log_scale(dataStruct(indShift+i_loop).chi_linear,...
                    dataStruct(indShift+i_loop).t_r, dt, st);
            else
                dataStruct(indShift+i_loop).chi_linear = ((16.0*sqrt(r_tip))/3)...
                    .* (dataStruct(indShift+i_loop).h_r_smooth.^1.5) ./ dataStruct(indShift+i_loop).Fdot;
                dataStruct(indShift+i_loop).chi_linear_log = log_scale(dataStruct(indShift+i_loop).chi_linear,...
                    dataStruct(indShift+i_loop).t_r_smooth, dt, st);
            end

            % Function Shape for Multiple Voigt with Linear Load Assumption
            [U_linear,F_linear_fit] = makeLinearViscoFunction(elasticSetting,fluidSetting,n_terms,dataStruct(indShift+i_loop).t_r,dt,st,fitLog);

            n_loss_fit = n_terms;

            % Make the label for this iteration
            switch avAnalysis
                case 'i'
                    fprintf(fid,'File No: %d\r\n',i_loop);
                case 'a'
                    fprintf(fid,'Approach Velocity: %4.3g nm/s\r\n',v_unique(i_loop));
            end

            if timeOperation
                timeStoreTemp = []; % Used to store the time for this loop
            end

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
                if ~advancedLog
                    F_log = log_scale(dataStruct(indShift+i_loop).F_r,dataStruct(indShift+i_loop).t_r,dt,st);
                    t_log = log_scale(dataStruct(indShift+i_loop).t_r,dataStruct(indShift+i_loop).t_r,dt,st);
                    h_log = log_scale(dataStruct(indShift+i_loop).h_r,dataStruct(indShift+i_loop).t_r,dt,st);
                else
                    [F_log, ind_set] = log_scale_advanced(dataStruct(indShift+i_loop).F_r,dataStruct(indShift+i_loop).t_r,dt,st);
                    t_log = dataStruct(indShift+i_loop).t_r(ind_set)';
                    h_log = dataStruct(indShift+i_loop).h_r(ind_set)';
                end

                F_schapery_area = dataStruct(indShift+i_loop).F_r;
                negativeIndexes = dataStruct(indShift+i_loop).F_r <=  0;
                x_fit = [dataStruct(indShift+i_loop).t_r, dataStruct(indShift+i_loop).F_r];
                x_fit(negativeIndexes,2) = 1e-25;
                timeVec = dataStruct(indShift+i_loop).t_r;
                maxwellInputs = [dataStruct(indShift+i_loop).t_r(dataStruct(indShift+i_loop).t_r>=0), (dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>=0).^1.5)];
                maxwellTimeVec = dataStruct(indShift+i_loop).t_r(dataStruct(indShift+i_loop).t_r>=0);
                x_fit_area = [maxwellInputs(:,1), sqrt(2*r_tip*dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>=0) - (dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>=0).^2))./r_tip, dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>=0)];
            else
                if ~advancedLog
                    F_log = log_scale(dataStruct(indShift+i_loop).F_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);
                    t_log = log_scale(dataStruct(indShift+i_loop).t_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);
                    h_log = log_scale(dataStruct(indShift+i_loop).h_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);
                else
                    [F_log, ind_set] = log_scale_advanced(dataStruct(indShift+i_loop).F_r_smooth,dataStruct(indShift+i_loop).t_r_smooth,dt,st);
                    t_log = dataStruct(indShift+i_loop).t_r_smooth(ind_set)';
                    h_log = dataStruct(indShift+i_loop).h_r_smooth(ind_set)';
                end

                F_schapery_area = dataStruct(indShift+i_loop).F_r_smooth;
                negativeIndexes = dataStruct(indShift+i_loop).F_r_smooth <=  0;
                x_fit = [dataStruct(indShift+i_loop).t_r_smooth, dataStruct(indShift+i_loop).F_r_smooth];
                x_fit(negativeIndexes,2) = 1e-25;
                timeVec = dataStruct(indShift+i_loop).t_r_smooth;
                maxwellInputs = [dataStruct(indShift+i_loop).t_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0), (dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0).^1.5)];
                maxwellTimeVec = dataStruct(indShift+i_loop).t_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0);
                x_fit_area = [maxwellInputs(:,1), sqrt(2*r_tip*dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0) - (dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0).^2))./r_tip, dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>=0)];
            end

            % Generate a frequency array in log scale
            de0 = 2.0.*pi.*(1./timeVec(end));
            maxi = 2.0.*pi.*(1./mean(diff(timeVec)));
            omega = log_tw(de0,maxi);

            dataStruct(indShift+i_loop).h_norm = h_norm;
            dataStruct(indShift+i_loop).F_norm = F_norm;

            dataStruct(indShift+i_loop).h_log = h_log;
            h_norm_log = log_scale(h_norm,timeVec,dt,st);
            dataStruct(indShift+i_loop).h_norm_log = h_norm_log;
            F_norm_log = log_scale(F_norm,timeVec,dt,st);
            dataStruct(indShift+i_loop).F_norm_log = F_norm_log;

            % For Schapery Maxwell
            dataStruct(indShift+i_loop).F_schapery_area = F_schapery_area;
            F_log_schapery_area = log_scale(F_schapery_area,timeVec,dt,st);
            dataStruct(indShift+i_loop).F_log_schapery_area = F_log_schapery_area;

            % Create the Generalized Voigt model
            [U_func,F_conv,F_conv_wrapper,lb,ub,subref,selector] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);

            % To recreate this data we will need these values later
            dataStruct(indShift+i_loop).nu_sample = nu_sample;
            dataStruct(indShift+i_loop).nu_tip = nu_tip;
            dataStruct(indShift+i_loop).E_tip = E_tip;
            dataStruct(indShift+i_loop).r_tip = r_tip;

            % Make Schapery Functions
            [~,E_star,A_conv,h_conv_schapery_area,h_conv_schapery_area_wrapper] = makeSchaperyModel(elasticSetting_maxwell,fluidSetting_maxwell,'n',maxwellTimeVec,dt,st,minTimescale,'maxwell',n_loss_fit,advancedLog,forwardFitTimescale,r_tip,E_tip,nu_tip,nu_sample);

            % Create the Gen. Maxwell Force Convolution Model
            [G_func_maxwell,Q_func_maxwell,h_conv_maxwell,h_conv_maxwell_wrapper,padSize] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_loss_fit,maxwellTimeVec,dt,st);

            % Create the Harmonic Gen. Maxwell Models
            [E_storage,E_loss,lossAngle,harmonic_wrapper,lb_maxwell,ub_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_loss_fit,maxwellTimeVec,minTimescale,padSize,forwardFitTimescale);

            % Open Parallel Pool of MATLAB Workers
            if isempty(gcp('nocreate'))
                %parpool(N_workers)
                parpool(8,'IdleTimeout', Inf)
            end

            if timeOperation
                fitSwitchStartTime(iVoigt,i_loop) = toc;
            end

            % Begin fitting
            switch lower(fitMethod)
                case 'iterative'
                    linearParams_loop = {};
                    maxwellParams_loop = {};
                    convParams_loop = {};
                    schaperyParams_loop = {};

                    if strcmp(elasticSetting,'y')
                        % Get the best elastic parameter for our iterative fit
                        [~,~,~,lb_elastic,ub_elastic,subref_elastic,selector_elastic] = makeGeneralizedVoigtModel('y','n',0,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);

                        x_fit_elastic = [log_scale(x_fit(:,1),x_fit(:,1),dt,st); log_scale(x_fit(:,2),x_fit(:,1),dt,st)]';
                        F_conv_wrapper_elastic = @(c,inputs) c(1).*inputs(:,2);

                        y_fit_elastic = (h_norm_log');
        %                     y_fit_elastic = (h_norm_log./(10.^floor(max(log10(h_norm_log))-min(log10(h_norm_log)))))';

                        if fitElasticFirst
                            beta_dist_elastic = zeros(length(ub_elastic),n_samples);
                            beta0_dist_elastic = zeros(size(beta_dist_elastic));
                            resnorm_dist_elastic = zeros(1,n_samples);

                            if timeOperation
                                elasticParallelStartTime(iVoigt,i_loop) = toc;
                                ticBytes(gcp);
                            end

                            ub_rand_elastic = log10(ub_elastic(1))-1;
        %                     ub_rand_elastic = -6;
                            lb_rand_elastic = log10(lb_elastic(1))+1;
                            beta0_elastic_array = logspace(ub_rand_elastic,lb_rand_elastic,n_samples);

                            parfor i = 1:n_samples
                                lsqoptions = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
                                    'MaxFunctionEvaluations', n_maxIterations,...
                                    'MaxIterations', n_maxIterations,...
                                    'FiniteDifferenceType','central',...
                                    'FunctionTolerance', 0,...
                                    'OptimalityTolerance', 0,...
                                    'StepTolerance', scaling,...
                                    'Display', 'none');

                                % Perform the fit using fmincon
        %                         Fsumsquares_elastic = @(c) sum(( abs((F_conv_wrapper_elastic(c,x_fit) - y_fit_elastic)) ./ abs(y_fit_elastic)).^2);
        %                         opts = optimoptions('fmincon','Algorithm','sqp',...
        %                             'MaxFunctionEvaluations', n_maxIterations,...
        %                             'MaxIterations', n_maxIterations,...
        %                             'Display','none',...
        %                             'FiniteDifferenceType','central',...
        %                             'FunctionTolerance', scaling,...
        %                             'OptimalityTolerance', scaling,...
        %                             'StepTolerance', scaling);
        % 
        %                         [betatemp,ressquared,eflag,outputu] = ...
        %                             fmincon(Fsumsquares_elastic,beta0_elastic_array(i),[],[],[],[],lb_elastic,ub_elastic,[],opts);


                                [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                    lsqcurvefit(F_conv_wrapper_elastic,beta0_elastic_array(i),x_fit_elastic,y_fit_elastic,lb_elastic,ub_elastic,lsqoptions);

                                % Test the fit (if necessary)
        %                         figure
        %                         scatter(x_fit_elastic(:,1),y_fit_elastic,'ro')
        %                         hold on
        %                         plot(x_fit_elastic(:,1),F_conv_wrapper_elastic(beta0_elastic_array(i),x_fit_elastic),'b-')
        %                         plot(x_fit_elastic(:,1),F_conv_wrapper_elastic(betatemp,x_fit_elastic),'r-')
        %                         set(gca,'XScale','log')
        %                         set(gca,'YScale','log')
        %                         hold off

                                beta_dist_elastic(:,i) = betatemp;
                                beta0_dist_elastic(:,i) = beta0_elastic_array(i);

    %                             resnorm_dist_elastic(i) = sum((residual./y_fit_loop).^2);
    %                             resnorm_dist_elastic(i) = sum((residual).^2);
    %                             resnorm_dist_elastic(i) = resnorm;
                                if forceResidual
                                    yModelIndentation = (1./((8*sqrt(r_tip))/(3*(1-nu_sample))).*F_conv_wrapper_elastic(betatemp,x_fit_elastic)).^(2/3);
                                    resnorm_dist_elastic(i) = sum(((yModelIndentation-h_log').^2)./movvar(h_log',3))./(length(h_log')-length(betatemp)); % Standard Error of Regression (S)
                                else
                                    resnorm_dist_elastic(i) = sum(((F_conv_wrapper_elastic(betatemp,x_fit_elastic)-y_fit_elastic).^2)./movvar(y_fit_elastic,3))./(length(y_fit_elastic)-length(betatemp)); % Standard Error of Regression (S)
                                end

                            end

                            if timeOperation
                                elasticParallelEndTime(iVoigt,i_loop) = toc;
                                elasticParallelBytes(iVoigt,i_loop) = getfield(tocBytes(gcp),{1});

                                fprintf('\nThe Elastic Parameter parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                                    (elasticParallelEndTime(iVoigt,i_loop)-elasticParallelStartTime(iVoigt,i_loop))/60, elasticParallelBytes(iVoigt,i_loop)*1e-6);
                            end

                            resnorm_dist_elastic_temp = resnorm_dist_elastic;
                            resnorm_dist_elastic_temp(resnorm_dist_elastic_temp == 0) = NaN;
                            [~,idx] = min(resnorm_dist_elastic_temp,[],'omitnan');
                            bestElasticTerm = beta_dist_elastic(:,idx);

                        else
                            if timeOperation
                                elasticParallelStartTime(iVoigt,i_loop) = toc;
                                ticBytes(gcp);
                            end
                            if timeOperation
                                elasticParallelEndTime(iVoigt,i_loop) = toc;
                                elasticParallelBytes(iVoigt,i_loop) = getfield(tocBytes(gcp),{1});
                            end

                        end
                    end

                    if timeOperation
                        iterativeFitStartTime(iVoigt,i_loop) = toc;
                    end

                    % Perform the iterative fitting
                    for i = 1:n_terms

                        [~,~,~,~,ub_loop,~,~] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,i,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);
                        [~,~,~,~,padSize_loop] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellTimeVec,dt,st);
                        [~,~,~,~,~,ub_loop_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellTimeVec,minTimescale,padSize_loop,forwardFitTimescale);

                        [tauInds,modulusInds] = getParamIndices(ub_loop,elasticSetting,fluidSetting);
                        [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(ub_loop_maxwell,elasticSetting_maxwell,fluidSetting_maxwell);
                        
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

                        if smoothData == 'n'
                            if clipData
                                if fitLog
                                    y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds))),timeVec,dt,st);
                                    y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec,dt,st);
                                    y_fit_linear_loop = log_scale(dataStruct(indShift+i_loop).chi_linear(x_fit(:,1) < max(ub_loop(tauInds))),timeVec,dt,st);
                                    y_fit_area_loop = log_scale(dataStruct(indShift+i_loop).F_schapery_area(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec,dt,st);
                                else
                                    y_fit_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                    y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                                    y_fit_linear_loop = dataStruct(indShift+i_loop).chi_linear(x_fit(:,1) < max(ub_loop(tauInds)));
                                    y_fit_area_loop = dataStruct(indShift+i_loop).F_schapery_area(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                                end
                                
                                maxwellInputs_loop = maxwellInputs(x_fit(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                x_fit_loop = x_fit(x_fit(:,1) < max(ub_loop(tauInds)),:);
                                x_fit_area_loop = x_fit_area(x_fit_area(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                                h_r_loop = dataStruct(indShift+i_loop).h_r(x_fit(:,1) < max(ub_loop(tauInds)));
                                h_norm_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                F_r_loop = dataStruct(indShift+i_loop).F_r(x_fit(:,1) < max(ub_loop(tauInds)));
                                F_norm_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                                F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                                F_log_loop = log_scale(F_r_loop,x_fit_area_loop(:,1),dt,st);
                            else
                                if fitLog
                                    y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm,timeVec,dt,st);
                                    y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm,timeVec,dt,st);
                                    y_fit_linear_loop = log_scale(dataStruct(indShift+i_loop).chi_linear,timeVec,dt,st);
                                    y_fit_area_loop = log_scale(dataStruct(indShift+i_loop).F_schapery_area,timeVec,dt,st);
                                else
                                    y_fit_loop = dataStruct(indShift+i_loop).h_norm;
                                    y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm;
                                    y_fit_linear_loop = dataStruct(indShift+i_loop).chi_linear;
                                    y_fit_area_loop = dataStruct(indShift+i_loop).F_schapery_area;
                                end
                                
                                maxwellInputs_loop = maxwellInputs;
                                x_fit_loop = x_fit;
                                x_fit_area_loop = x_fit_area;
                                t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                                h_r_loop = dataStruct(indShift+i_loop).h_r;
                                h_norm_loop = dataStruct(indShift+i_loop).h_norm;
                                F_r_loop = dataStruct(indShift+i_loop).F_r;
                                F_norm_loop = dataStruct(indShift+i_loop).F_norm;
                                h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                                F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                                F_log_loop = log_scale(F_r_loop,x_fit_area_loop(:,1),dt,st);
                            end
                        else
                            if clipData
                                if fitLog
                                    y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds))),timeVec,dt,st);
                                    y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec,dt,st);
                                    y_fit_linear_loop = log_scale(dataStruct(indShift+i_loop).chi_linear(x_fit(:,1) < max(ub_loop(tauInds))),timeVec,dt,st);
                                    y_fit_area_loop = log_scale(dataStruct(indShift+i_loop).F_schapery_area(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec,dt,st);
                                else
                                    y_fit_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                    y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                                    y_fit_linear_loop = dataStruct(indShift+i_loop).chi_linear(x_fit(:,1) < max(ub_loop(tauInds)));
                                    y_fit_area_loop = dataStruct(indShift+i_loop).F_schapery_area(x_fit(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                                end
                                
                                maxwellInputs_loop = maxwellInputs(x_fit(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                x_fit_loop = x_fit(x_fit(:,1) < max(ub_loop(tauInds)),:);
                                x_fit_area_loop = x_fit_area(x_fit_area(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                                h_r_loop = dataStruct(indShift+i_loop).h_r_smooth(x_fit(:,1) < max(ub_loop(tauInds)));
                                h_norm_loop = dataStruct(indShift+i_loop).h_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                F_r_loop = dataStruct(indShift+i_loop).F_r_smooth(x_fit(:,1) < max(ub_loop(tauInds)));
                                F_norm_loop = dataStruct(indShift+i_loop).F_norm(x_fit(:,1) < max(ub_loop(tauInds)));
                                h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                                F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                                F_log_loop = log_scale(F_r_loop,x_fit_area_loop(:,1),dt,st);
                            else
                                if fitLog
                                    y_fit_loop = log_scale(dataStruct(indShift+i_loop).h_norm,timeVec,dt,st);
                                    y_fit_maxwell_loop = log_scale(dataStruct(indShift+i_loop).F_norm,timeVec,dt,st);
                                    y_fit_linear_loop = log_scale(dataStruct(indShift+i_loop).chi_linear,timeVec,dt,st);
                                    y_fit_area_loop = log_scale(dataStruct(indShift+i_loop).F_schapery_area,timeVec,dt,st);
                                else
                                    y_fit_loop = dataStruct(indShift+i_loop).h_norm;
                                    y_fit_maxwell_loop = dataStruct(indShift+i_loop).F_norm;
                                    y_fit_linear_loop = dataStruct(indShift+i_loop).chi_linear;
                                    y_fit_area_loop = dataStruct(indShift+i_loop).F_schapery_area;
                                end
                                
                                maxwellInputs_loop = maxwellInputs;
                                x_fit_loop = x_fit;
                                x_fit_area_loop = x_fit_area;
                                t_log_loop = log_scale(x_fit_loop(:,1),timeVec,dt,st);
                                h_r_loop = dataStruct(indShift+i_loop).h_r_smooth;
                                h_norm_loop = dataStruct(indShift+i_loop).h_norm;
                                F_r_loop = dataStruct(indShift+i_loop).F_r_smooth;
                                F_norm_loop = dataStruct(indShift+i_loop).F_norm;
                                h_norm_log_loop = log_scale(h_norm_loop,x_fit_loop(:,1),dt,st);
                                F_norm_log_loop = log_scale(F_norm_loop,x_fit_loop(:,1),dt,st);
                                F_log_loop = log_scale(F_r_loop,x_fit_area_loop(:,1),dt,st);
                            end
                        end
                        
                        if fitLog
                            t_plot_loop = t_log_loop;
                        else
                            t_plot_loop = x_fit_loop(:,1);
                        end

                        h_log_loop = log_scale(h_r_loop,x_fit_loop(:,1),dt,st);
                        F_schapery_area_loop = F_r_loop;

                        % Make the Generalized Voigt Model for this loop
                        [U_func_loop,F_conv_loop,F_conv_wrapper_loop,lb_loop,ub_loop,subref_loop,selector_loop] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,i,x_fit_loop(:,1),dt,st,minTimescale,advancedLog,forwardFitTimescale);

                        % Make the Schapery functions for this loop
                        [Q_func_loop,E_star_loop,A_conv_loop,h_conv_schapery_area_loop,h_conv_schapery_area_wrapper_loop] = makeSchaperyModel(elasticSetting_maxwell,fluidSetting_maxwell,'n',maxwellInputs_loop(:,1),dt,st,minTimescale,'maxwell',i,advancedLog,forwardFitTimescale,r_tip,E_tip,nu_tip,nu_sample);

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing the Area vs L&R
%                         figure
%                         plot(timeVec,((8*sqrt(r_tip))/(3*(1-nu_sample))).*ones(size(timeVec)))
%                         hold on
%                         plot(x_fit_area(:,1),A_conv_loop(0,x_fit_area))
%                         legend({'L&R','Schapery Area Norm'},'location','southoutside')
%                         limSet = findobj(gca, '-property', 'ydata');
%                         limSet = get(limSet, 'YData');
%                         limSet = [limSet{:}];
%                         set(gca,'yscale','log')
%                         set(gca,'xscale','log')
%                         hold off
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Testing the Norm Data
%                         figure
%                         scatter(timeVec,F_schapery_area)
%                         hold on
%                         plot(timeVec,h_conv_schapery_area_loop(manualParams,x_fit_area_loop))
%                         legend({'Data','Model'},'location','southoutside')
%                         limSet = findobj(gca, '-property', 'ydata');
%                         limSet = get(limSet, 'YData');
%                         limSet = [limSet{:}];
%                         set(gca,'ylim',[min(limSet,[],'all')/limMargin max(limSet,[],'all')*limMargin])
%                         set(gca,'yscale','log')
%                         set(gca,'xscale','log')
%                         hold off
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % Function Shape for Multiple Voigt with Linear Load Assumption
                        [U_linear_loop,F_linear_fit_loop] = makeLinearViscoFunction(elasticSetting,fluidSetting,i,x_fit_loop(:,1),dt,st,fitLog);

                        % Create the Gen. Maxwell Force Convolution Model
                        [G_func_maxwell_loop,Q_func_maxwell_loop,h_conv_maxwell_loop,h_conv_maxwell_wrapper_loop,padSize_loop] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellInputs_loop(:,1),dt,st);

                        % Create the Harmonic Gen. Maxwell Models
                        [E_storage_loop,E_loss_loop,lossAngle_loop,harmonic_wrapper_loop,lb_loop_maxwell,ub_loop_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellInputs_loop(:,1),minTimescale,padSize_loop,forwardFitTimescale);

                        beta0_linear = NaN(size(ub_loop));
                        beta0_maxwell = NaN(size(ub_loop_maxwell));
                        beta0_conv = NaN(size(ub_loop));
                        beta0_schapery = NaN(size(ub_loop_maxwell));

                        ub_randLimit_multigrid = 1*ones(size(ub_loop));
                        ub_randLimit_multigrid(tauInds) = ub_loop(tauInds);

                        lb_randLimit_multigrid = -5*ones(size(lb_loop));
                        lb_randLimit_multigrid(tauInds) = lb_loop(tauInds);

                        ub_randLimit_maxwell_multigrid = 5*ones(size(ub_loop_maxwell));
                        ub_randLimit_maxwell_multigrid(tauInds_maxwell) = ub_loop_maxwell(tauInds_maxwell);

                        lb_randLimit_maxwell_multigrid = -1*ones(size(lb_loop_maxwell));
                        lb_randLimit_maxwell_multigrid(tauInds_maxwell) = lb_loop_maxwell(tauInds_maxwell);

                        ub_randLimit_schapery_multigrid = ub_randLimit_maxwell_multigrid;
                        lb_randLimit_schapery_multigrid = lb_randLimit_maxwell_multigrid;

                        % For the Voigt models
                        if i > 1
                            if fluidSetting == 'y'
                                beta0_linear(1:length(linearParams_loop{i-1})-1) = linearParams_loop{i-1}(1:end-1);
                                beta0_conv(1:length(convParams_loop{i-1})-1) = convParams_loop{i-1}(1:end-1);

                                beta0_linear(end) = 0;
                                beta0_conv(end) = 0;
                            else
                                beta0_linear(1:length(linearParams_loop{i-1})) = linearParams_loop{i-1};
                                if strcmp(elasticSetting,'y') && isnan(beta0_linear(1))
                                    beta0_linear(1) = bestElasticTerm;
                                end
                                beta0_conv(1:length(convParams_loop{i-1})) = convParams_loop{i-1};
                                if strcmp(elasticSetting,'y') && isnan(beta0_conv(1))
                                    beta0_conv(1) = bestElasticTerm;
                                end
                            end
                        else
                            % Manual override beta0 for random sampling
                            beta0_loop = ((ub_loop-lb_loop).*rand(size(ub_loop)) + lb_loop);
                            beta0_linear = beta0_loop;
                            beta0_conv = beta0_loop;
                        end

                        % For the Maxwell models
                        if i > 1
                            if fluidSetting_maxwell == 'y'
                                beta0_maxwell(1:length(maxwellParams_loop{i-1})-1) = maxwellParams_loop{i-1}(1:end-1);
                                beta0_schapery(1:length(schaperyParams_loop{i-1})-1) = schaperyParams_loop{i-1}(1:end-1);

                                beta0_maxwell(end) = 0;
                                beta0_schapery(end) = 0;
                            else
                                beta0_maxwell(1:length(maxwellParams_loop{i-1})) = maxwellParams_loop{i-1};
                                if strcmp(elasticSetting_maxwell,'y') && isnan(beta0_maxwell(1))
                                    beta0_maxwell(1) = 1/bestElasticTerm;
                                end
                                beta0_schapery(1:length(schaperyParams_loop{i-1})) = schaperyParams_loop{i-1};
                                if strcmp(elasticSetting_maxwell,'y') && isnan(beta0_schapery(1))
                                    beta0_schapery(1) = 1/bestElasticTerm;
                                end
                            end
                        else
                            % Manual override beta0 for random sampling
                            beta0_loop = ((ub_loop_maxwell-lb_loop_maxwell).*rand(size(ub_loop_maxwell)) + lb_loop_maxwell);
                            beta0_maxwell = beta0_loop;
                            beta0_schapery = beta0_loop;
                        end

                        betaGrid = [];
                        betaGrid_maxwell = [];
                        betaGrid_schapery = [];
                        betaGridIn = {};
                        betaGridIn_maxwell = {};
                        betaGridIn_schapery = {};

                        if bruteForceBeta
                            % Create N-Dimensional Grid of guesses
                            tauStep = 2; % Number of steps in the range of tau starting values (set as low as possible).
                            newInds_conv = isnan(beta0_conv);
                            newInds_maxwell = isnan(beta0_maxwell);
                            newInds_schapery = isnan(beta0_schapery);

                            [tauInds,modulusInds] = getParamIndices(ub_loop,elasticSetting,fluidSetting);
                            [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(ub_loop_maxwell,elasticSetting_maxwell,fluidSetting_maxwell);
                            
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

                            numStepsPerTerm = floor((n_samples-length(tauInds)*tauStep)^(1./(length(ub_loop(modulusInds))-elasticInd)))/(length(tauInds)*tauStep);

                            for kk = 1:(length(ub_loop))
                                if any(kk == find(newInds_conv,1)) || isempty(find(newInds_conv,1))
                                    if (strcmp(elasticSetting,'n') || kk > elasticInd) && fitElasticFirst
                                        if any(kk == tauInds)
                                            loopStep = tauStep;
                                        elseif any(kk == modulusInds)
                                            loopStep = numStepsPerTerm;
                                        end
                                        betaGridIn{kk} = logspace(log10(lb_loop(kk)),log10(ub_loop(kk)),loopStep);
                                    else
                                        betaGridIn{kk} = bestElasticTerm;
                                    end
                                else
                                    betaGridIn{kk} = convParams_loop{i-1};
                                end
                            end

                            for kk = 1:(length(ub_loop_maxwell))
                                if any(kk == find(newInds_maxwell,1)) || isempty(find(newInds_maxwell,1))
                                    if (strcmp(elasticSetting_maxwell,'n') || kk > elasticInd) && fitElasticFirst
                                        if any(kk == tauInds)
                                            loopStep = tauStep;
                                        elseif any(kk == modulusInds)
                                            loopStep = numStepsPerTerm;
                                        end
                                        betaGridIn_maxwell{kk} = logspace(log10(lb_loop_maxwell(kk)),log10(ub_loop_maxwell(kk)),loopStep);
                                        betaGridIn_schapery{kk} = logspace(log10(lb_loop_maxwell(kk)),log10(ub_loop_maxwell(kk)),loopStep);
                                    else
                                        betaGridIn_maxwell{kk} = 1./bestElasticTerm;
                                        betaGridIn_schapery{kk} = 1./bestElasticTerm;
                                    end
                                else
                                    betaGridIn_maxwell{kk} = maxwellParams_loop{i-1};
                                    betaGridIn_schapery{kk} = schaperyParams_loop{i-1};
                                end
                            end

                            c = betaGridIn;
                            [c{:}] = ndgrid(c{:});
                            nspecs = length(c);
                            betaGrid = reshape(cat(nspecs+1,c{:}),[],nspecs);

                            c = betaGridIn_maxwell;
                            [c{:}] = ndgrid(c{:});
                            nspecs = length(c);
                            betaGrid_maxwell = reshape(cat(nspecs+1,c{:}),[],nspecs);

                            c = betaGridIn_schapery;
                            [c{:}] = ndgrid(c{:});
                            nspecs = length(c);
                            betaGrid_schapery = reshape(cat(nspecs+1,c{:}),[],nspecs);

                            loopLim = min([size(betaGrid,1), size(betaGrid_maxwell,1), size(betaGrid_schapery,1)]);
                            fprintf('\nBrute Force Grid Size: Iterations Set to %d\n',loopLim);

                        else

                            % Need to create placeholders that look real
                            % because otherwise the parfor pre-processing will
                            % fail.
                            loopLim = n_samples;
                            betaGrid = ones(loopLim,length(ub_loop));
                            betaGrid_maxwell = ones(loopLim,length(ub_loop_maxwell));
                            betaGrid_schapery = ones(loopLim,length(ub_loop_maxwell));

                        end

                        beta_dist = zeros(length(ub_loop),loopLim);
                        beta0_dist = zeros(size(beta_dist));
                        resnorm_dist = zeros(1,loopLim);
                        beta_dist_maxwell = zeros(length(ub_loop_maxwell),loopLim);
                        beta0_dist_maxwell = zeros(size(beta_dist_maxwell));
                        resnorm_dist_maxwell = zeros(1,loopLim);
                        beta_dist_schapery_area = zeros(length(ub_loop_maxwell),loopLim);
                        beta0_dist_schapery_area = zeros(size(beta_dist_schapery_area));
                        resnorm_dist_schapery_area = zeros(1,loopLim);
                        beta_dist_linear = zeros(length(ub_loop),loopLim);
                        beta0_dist_linear = zeros(size(beta_dist_linear));
                        resnorm_dist_linear = zeros(1,loopLim);

                        if timeOperation
                            parallelFitStartTime(iVoigt,i_loop,i) = toc;
                            convTime = Par(loopLim);
                            maxwellTime = Par(loopLim);
                            schaperyTime = Par(loopLim);
                            if i == 1
                               convTimeTotal = [];
                               maxwellTimeTotal = [];
                               schaperyTimeTotal = [];
                            end
                            ticBytes(gcp);
                            if ~exist('userMem','var')
                                %[userMem,systemMem] = memory;
                                userMem = 0;
                                systemMem = 0;
                            end
                            ramused = NaN(1,loopLim);
                            ramav = NaN(1,loopLim);
                            ramtime = NaN(1,loopLim);
                        end

                        errorCalcVoigt = @(c,inputs) sum(((F_conv_wrapper_loop(c,inputs)-y_fit_loop).^2)./movvar(y_fit_loop,3))./(length(y_fit_loop)-length(c)); % Standard Error of Regression (S)
                        errorCalcMaxwell = @(c,inputs) sum(((h_conv_maxwell_wrapper_loop(c,inputs)-y_fit_maxwell_loop).^2)./movvar(y_fit_maxwell_loop,3))./(length(y_fit_maxwell_loop)-length(c)); % Standard Error of Regression (S)
                        errorCalcSchapery = @(c,inputs) sum(((h_conv_schapery_area_wrapper_loop(c,inputs)-y_fit_area_loop).^2)./movvar(y_fit_area_loop,3))./(length(y_fit_area_loop)-length(c)); % Standard Error of Regression (S)

                        switch lower(optimizationMethod)
                            case 'gd'

                                printData = 0;
                                plotCost = 0;
                                plotFit = 0;
                                n_plot = 1000;

                                integralCalc = @(c,inputs,dQdI) subref_loop(convnfft(dQdI(c,inputs), inputs(:,2),'full'))*dt;
                                jacobianCalc = @(c,inputs,dQdI) (dQdI(c,inputs));

                                if plotCost
                                    if exist('costPlot','var')
                                        clearvars costPlot
                                    end
                                    if exist('costPlot','var')
                                        figure(costPlot);
                                        clf;
                                    else
                                        costPlot = figure('position',[25 25 600 400]);
                                    end
                                    if exist('fitPlot','var')
                                        clearvars fitPlot
                                    end
                                    if exist('fitPlot','var')
                                        figure(fitPlot);
                                        clf;
                                    else
                                        fitPlot = figure('position',[125 25 600 400]);
                                    end
                                end

                                ydata_loop = F_r_loop;
                                ydata_voigt_loop = h_r_loop;
                                inputs_loop = [x_fit_loop(:,1),...
                                        F_r_loop];
                                inputsForce = [x_fit_loop(:,1),...
                                    h_r_loop.^1.5];
                                inputsForce_schapery = [x_fit_loop(x_fit_loop(:,1)>0,1),...
                                    sqrt(2*r_tip*h_r_loop(x_fit_loop(:,1)>0)...
                                    - (h_r_loop(x_fit_loop(:,1)>0).^2))./r_tip];
                                timeVec_loop = x_fit_loop(:,1);
                                alphaModel = (8*sqrt(r_tip))/(3*(1-nu_sample));
                                if ~strcmp(GDmethod,'alrm-norm')
                                    y_fit_GD_loop = h_norm_loop;
                                    y_fit_maxwell_GD_loop = F_norm_loop;
                                    y_fit_schapery_GD_loop = F_schapery_area_loop;
                                else
                                    y_fit_GD_loop = subref_loop(convnfft(h_norm_loop,1./F_r_loop,'full'))*dt;
                                    y_fit_maxwell_GD_loop = subref_loop(convnfft(F_norm_loop,1./(h_r_loop.^1.5),'full'))*dt;
                                    y_fit_schapery_GD_loop = subref_loop(convnfft(F_schapery_area_loop, ...
                                        1./((sqrt(2*r_tip*h_r_loop(x_fit_loop(:,1)>0) -...
                                        (h_r_loop(x_fit_loop(:,1)>0).^2))./r_tip) .*...
                                        (2*pi*r_tip*h_r_loop(x_fit_loop(:,1)>0))),'full'))*dt;
                                end

                                if symbolicJacobian
                                    inputs_loop = horzcat(inputs_loop,y_fit_GD_loop);
                                    inputsForce = horzcat(inputsForce,y_fit_maxwell_GD_loop);
                                    inputsForce_schapery = horzcat(inputsForce_schapery,y_fit_schapery_GD_loop);
                                end

                                multiGridFlag = 1;
                                multiGrid_iters = 1;

                                while multiGridFlag

                                    parfor j = 1:loopLim
                                        warning('off');
                                        % Initialize
                                        newInds_linear = [];
                                        newInds_conv = [];
                                        newInds_maxwell = [];
                                        newInds_schapery = [];
                                        beta0_linear_loop = [];
                                        beta0_conv_loop = [];
                                        beta0_maxwell_loop = [];
                                        beta0_schapery_loop = [];
                                        ub_randLimit = [];
                                        lb_randLimit = [];
                                        ub_randLimit_maxwell = [];
                                        lb_randLimit_maxwell = [];
                                        ub_randLimit_schapery = [];
                                        lb_randLimit_schapery = [];
                                        beta0_temp_loop = [];
                                        limScale = [];
                                        tauInds = [];
                                        oldTau = [];
                                        modulusInds = [];
                                        oldModuli = [];
                                        newParams_temp = [];

                                        % Find the new indices
                                        newInds_linear = isnan(beta0_linear);
                                        newInds_conv = isnan(beta0_conv);
                                        newInds_maxwell = isnan(beta0_maxwell);
                                        newInds_schapery = isnan(beta0_schapery);

                                        % Grab the base values for the beta0 cases
                                        beta0_linear_loop = beta0_linear;
                                        beta0_conv_loop = beta0_conv;
                                        beta0_maxwell_loop = beta0_maxwell;
                                        beta0_schapery_loop = beta0_schapery;

                                        % Control the bounds AND the initial parameters
                                        ub_randLimit = ub_randLimit_multigrid;
                                        lb_randLimit = lb_randLimit_multigrid;
                                        ub_randLimit_maxwell = ub_randLimit_maxwell_multigrid;
                                        lb_randLimit_maxwell = lb_randLimit_maxwell_multigrid;
                                        ub_randLimit_schapery = ub_randLimit_schapery_multigrid;
                                        lb_randLimit_schapery = lb_randLimit_schapery_multigrid;
                                        ub_fluidityRandLimit_init = 1e-20;
                                        limScale = getfield(logspace(5,0,n_terms),{i})*relaxationFactor;

                                        % Include or remove elastic term from bound
                                        % loosening that happens below. Use 1 for
                                        % inclusion, empty brackets [] for ignoring.
                                        elasticInd = 1;

            %                             oldModReductionFactor = limScale*0.5;
                                        oldModReductionFactor = 1;
            %                             oldModReductionFactor = 1/relaxationFactor;

                                        errorarray = NaN(1,n_maxIterations);

                                        switch GDmethod
                                            case {'standard','adadelta','sgd'}
                                                alphaarray = NaN(length(beta0_conv),n_maxIterations);
                                            case {'alrm','alrm-norm'}
                                                alphaarray = NaN(1,n_maxIterations);
                                        end

                                        % Create random starting point
                                        [beta0_conv_loop,tauInds,modulusInds] = makeRandomParams(beta0_conv,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,newInds_conv);
                                        
                                        weightArray = ones(size(beta0_conv_loop));
                                        if scaleParams
                                            numTau = length(tauInds);
                                            if strcmp(fluidSetting,'y')
                                               numTau = numTau-1;
                                            end
                                            if i == 1
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            F_conv_loop_weight = logscaleParams(F_conv_loop);
                                                        else
                                                            F_conv_loop_weight = logscaleParams(U_func_loop);
                                                        end

                                                    otherwise
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(modulusInds) = 10.^((max(ub_randLimit(modulusInds))-min(lb_randLimit(modulusInds)))/2);
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            F_conv_loop_weight = weightParams(F_conv_loop,weightArray);
                                                        else
                                                            F_conv_loop_weight = weightParams(U_func_loop,weightArray);
                                                        end
                                                end
                                            else
                                                switch lower(paramScale)
                                                    case 'equal'
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_conv(intersect(find(~newInds_conv),modulusInds))))));
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            F_conv_loop_weight = weightParams(F_conv_loop,weightArray);
                                                        else
                                                            F_conv_loop_weight = weightParams(U_func_loop,weightArray);
                                                        end
                                                    case 'individual'
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(intersect(find(~newInds_conv),modulusInds)) = 10.^(floor(log10(beta0_conv(intersect(find(~newInds_conv),modulusInds)))));
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            F_conv_loop_weight = weightParams(F_conv_loop,weightArray);
                                                        else
                                                            F_conv_loop_weight = weightParams(U_func_loop,weightArray);
                                                        end
                                                    case 'logscale'
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            F_conv_loop_weight = logscaleParams(F_conv_loop);
                                                        else
                                                            F_conv_loop_weight = logscaleParams(U_func_loop);
                                                        end
                                                end
                                            end
                                        else
                                            if ~strcmp(GDmethod,'alrm-norm')
                                                F_conv_loop_weight = weightParams(F_conv_loop,weightArray);
                                            else
                                                F_conv_loop_weight = weightParams(U_func_loop,weightArray);
                                            end
                                        end

                                        if plotCost
                                            legendEntries = {};
                                            for jj = 1:length(beta0_conv_loop)
                                                if ismember(jj,tauInds)
                                                    if strcmp(fluidSetting,'y') && jj == length(beta0_conv_loop)
                                                        legendEntries = horzcat(legendEntries,sprintf('$\\phi_f$'));
                                                    else
                                                        legendEntries = horzcat(legendEntries,sprintf('$\\tau_%d$',find(tauInds==jj)));
                                                    end
                                                elseif ismember(jj,modulusInds)
                                                    if strcmp(elasticSetting,'y') && jj == 1
                                                        legendEntries = horzcat(legendEntries,sprintf('$J_g$'));
                                                    else
                                                        legendEntries = horzcat(legendEntries,sprintf('$J_%d$',find(modulusInds==jj)));
                                                    end
                                                end
                                            end
                                        end

                                        iters = 0;
                                        endFlag = 0;
                                        perturbWait = 0;

                                        % Gradient Descent
                                        if timeOperation
                                            Par.tic;
                                        end

                                        while ~endFlag

                                            iters = iters + 1;
                                            if iters > 1
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    currentParams = newParams_voigt./weightArray;
                                                else
                                                    currentParams = log10(newParams_voigt);
                                                end
                                            end
                                            currentDiff = zeros(size(beta0_conv_loop));
                                            if iters == 1
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    currentParams = beta0_conv_loop./weightArray;
                                                else
                                                    currentParams = log10(beta0_conv_loop);
                                                end
                                                alphaTemp = alphaV;
                                                oldStep = zeros(size(currentDiff));
                                                gammaTemp = gammaV;
                                                runningAvGradient = zeros(size(currentDiff));
                                                currentUpdate = zeros(size(currentDiff));
                                                runningAvUpdate = zeros(size(currentDiff));
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    initialError = errorCalcVoigt(currentParams.*weightArray,x_fit_loop); % Standard Error of Regression (S)
                                                else
                                                    initialError = errorCalcVoigt(10.^currentParams,x_fit_loop); % Standard Error of Regression (S)
                                                end
                                                if ~symbolicJacobian
                                                    dQdI = cell(size(beta0_conv_loop));
                                                    for jj = 1:length(beta0_conv_loop)
                                                        if jj == 1 && strcmp(elasticSetting,'y')
        %                                                     dQdI{jj} = @(c,inputs) ones(size(inputs(:,1)));
                                                            dQdI{jj} = @(c,inputs) [1;zeros(length(inputs(:,1))-1,1)];
                                                        elseif any(jj == modulusInds)
                                                            tempFunc = sprintf('@(c,inputs) (1./c(%d)) .* exp( (-inputs(:,1))./c(%d) )',jj+1,jj+1);
                                                            dQdI{jj} = str2func(tempFunc);
                                                        elseif any(jj == tauInds)
                                                            tempFunc = sprintf('@(c,inputs) (c(%d)./(c(%d).^2)).*exp( (-inputs(:,1))./c(%d) ).*( (inputs(:,1)./c(%d)) - 1 )',jj-1,jj,jj,jj);
                                                            dQdI{jj} = str2func(tempFunc);
                                                        else
                                                            error('Error creating the derivative for term %d!',jj)
                                                        end
                                                    end
                                                else
                                                    [J_sym,H_sym,dQdI] = getObjectiveDerivatives(beta0_conv_loop,'gkv',elasticSetting,fluidSetting);
                                                end
                                            end

                                            % Calculate gradients inside integral term
                                            for jj = 1:length(beta0_conv_loop)
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        if ~normalizeGD
                                                            switch GDmethod
                                                                case 'standard'
                        %                                             currentDiff(jj) = (2/length(y_fit_GD_loop)).*sum((F_conv_loop(currentParams,x_fit_loop)-y_fit_GD_loop).*integralCalc(currentParams,x_fit_loop,dQdI{jj})); % Derivative of MSE
                                                                    currentDiff(jj) = (1/length(y_fit_GD_loop)).*sum(integralCalc(10.^currentParams,inputs_loop,dQdI{jj})-y_fit_GD_loop); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = (2/length(y_fit_GD_loop)).*sum((F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*integralCalc(10.^currentParams,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = (2/length(y_fit_GD_loop)).*sum((F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*jacobianCalc(10.^currentParams,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        else
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = sum((1+100.*inputs_loop(:,1).^2).*(integralCalc(10.^currentParams,inputs_loop,dQdI{jj})-y_fit_GD_loop)); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputs_loop(:,1).^2).*(F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*integralCalc(10.^currentParams,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputs_loop(:,1).^2).*(F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*jacobianCalc(10.^currentParams,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        end

                                                    otherwise
                                                        if ~normalizeGD
                                                            switch GDmethod
                                                                case 'standard'
                        %                                             currentDiff(jj) = (2/length(y_fit_GD_loop)).*sum((F_conv_loop(currentParams,x_fit_loop)-y_fit_GD_loop).*integralCalc(currentParams,x_fit_loop,dQdI{jj})); % Derivative of MSE
                                                                    currentDiff(jj) = (1/length(y_fit_GD_loop)).*sum(integralCalc(currentParams.*weightArray,inputs_loop,dQdI{jj})-y_fit_GD_loop); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = (2/length(y_fit_GD_loop)).*sum((F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*integralCalc(currentParams.*weightArray,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = (2/length(y_fit_GD_loop)).*sum((F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*jacobianCalc(currentParams.*weightArray,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        else
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = sum((1+100.*inputs_loop(:,1).^2).*(integralCalc(currentParams.*weightArray,inputs_loop,dQdI{jj})-y_fit_GD_loop)); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputs_loop(:,1).^2).*(F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*integralCalc(currentParams.*weightArray,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputs_loop(:,1).^2).*(F_conv_loop_weight(currentParams,inputs_loop)-y_fit_GD_loop).*jacobianCalc(currentParams.*weightArray,inputs_loop,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        end

                                                end
                                            end

                                            % Update current parameters
                                            switch GDmethod
                                                case 'standard'
                                                    newParams_voigt = currentParams - alphaTemp.*currentDiff;
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_voigt = 10.^newParams_voigt;
                                                        otherwise
                                                            newParams_voigt = newParams_voigt.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            newParams_voigt = currentParams + rand(size(currentParams)).*currentParams;
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop - newParams_voigt) < 0) ) || ( any((newParams_voigt - lb_loop) < 0) )
                                                        newParams_voigt((ub_loop - newParams_voigt) < 0) = ub_loop((ub_loop - newParams_voigt) < 0);
                                                        newParams_voigt((newParams_voigt - lb_loop) < 0) = lb_loop((newParams_voigt - lb_loop) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)
                                                    alphaarray(:,iters) = -alphaTemp.*currentDiff;

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Parameter Update vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Parameter Update')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case 'sgd'
                                                    alphaarray(:,iters) = -alphaTemp.*(currentDiff.*(sqrt(3)*(2*rand(size(currentDiff))-1)));
                                                    newParams_voigt = currentParams + alphaarray(:,iters)';
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_voigt = 10.^newParams_voigt;
                                                        otherwise
                                                            newParams_voigt = newParams_voigt.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            newParams_voigt = currentParams + rand(size(currentParams)).*currentParams;
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop - newParams_voigt) < 0) ) || ( any((newParams_voigt - lb_loop) < 0) )
                                                        newParams_voigt((ub_loop - newParams_voigt) < 0) = ub_loop((ub_loop - newParams_voigt) < 0);
                                                        newParams_voigt((newParams_voigt - lb_loop) < 0) = lb_loop((newParams_voigt - lb_loop) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Parameter Update vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Parameter Update')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case {'alrm','alrm-norm'}
                                                    newParams_voigt = currentParams + (gammaTemp.*oldStep - (1-gammaTemp).*alphaTemp.*currentDiff);
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_voigt = 10.^newParams_voigt;
                                                        otherwise
                                                            newParams_voigt = newParams_voigt.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            switch perturbMethod
                                                                case 'params'
                                                                    newParams_voigt = currentParams + rand(size(currentParams)).*currentParams;
                                                                case 'LR'
                                                                    alphaTemp = alphaV;
                                                            end
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop - newParams_voigt) < 0) ) || ( any((newParams_voigt - lb_loop) < 0) )
                                                        newParams_voigt((ub_loop - newParams_voigt) < 0) = ub_loop((ub_loop - newParams_voigt) < 0);
                                                        newParams_voigt((newParams_voigt - lb_loop) < 0) = lb_loop((newParams_voigt - lb_loop) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)

                                                    % Test if a valid step was taken
                                                    if iters == 1
                                                        errorChange = errorarray(iters) / initialError; % Error Ratio
                                                    else
                                                        errorChange = errorarray(iters) / errorarray(iters-1); % Error Ratio
                                                    end

                                                    if errorChange < zeta
                                                        alphaTemp = alphaTemp * eta; % Increase learning rate
                                                        gammaTemp = gammaV;
                                                        oldStep = (gammaTemp.*oldStep - (1-gammaTemp).*alphaTemp.*currentDiff);
                                                    else
                                                        alphaTemp = alphaTemp * rho; % Decrease learning rate
                                                        newParams_voigt = currentParams; % Throw away parameter updates
                                                        switch lower(paramScale)
                                                            case 'logscale'
                                                                newParams_voigt = 10.^newParams_voigt;
                                                            otherwise
                                                                newParams_voigt = newParams_voigt.*weightArray;
                                                        end
                                                        errorarray(iters) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)
                                                        gammaTemp = 0;
                                                    end
                                                    alphaarray(iters) = alphaTemp;

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        scatter((1:iters),alphaarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,alphaV,'rx')
                                                        title('Learning Parameter vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Learning Parameter')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case 'adadelta'
                                                    runningAvGradient = rho_adaV.*(runningAvGradient)+(1-rho_adaV).*(currentDiff.^2);
                                                    currentUpdate = -((sqrt(runningAvUpdate+epsilon_adadelta))./sqrt(runningAvGradient+epsilon_adadelta)).*currentDiff;
                                                    runningAvUpdate = rho_adaV.*(runningAvUpdate)+(1-rho_adaV).*(currentUpdate.^2);
                                                    newParams_voigt = currentParams + currentUpdate;
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_voigt = 10.^newParams_voigt;
                                                        otherwise
                                                            newParams_voigt = newParams_voigt.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            switch perturbMethod
                                                                case 'params'
                                                                    newParams_voigt = currentParams + rand(size(currentParams)).*currentParams;
                                                                case 'LR'
                                                                    runningAvUpdate = (runningAvUpdate./abs(runningAvUpdate)).*abs(runningAvUpdate).^(rand(size(runningAvUpdate))+1);
                                                            end
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop - newParams_voigt) < 0) ) || ( any((newParams_voigt - lb_loop) < 0) )
                                                        newParams_voigt((ub_loop - newParams_voigt) < 0) = ub_loop((ub_loop - newParams_voigt) < 0);
                                                        newParams_voigt((newParams_voigt - lb_loop) < 0) = lb_loop((newParams_voigt - lb_loop) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)
                                                    alphaarray(:,iters) = (sqrt(runningAvUpdate)./sqrt(runningAvGradient+epsilon_adadelta));

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Learning Parameter vs. Iteration No. [Voigt]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Learning Parameter')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                            end

                                            perturbWait = perturbWait + 1;

                                            if iters > 1
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        temp = 10.^currentParams;
                                                    otherwise
                                                        temp = currentParams.*weightArray;
                                                end
                                                if max(abs((temp-newParams_voigt)./temp)) < accuracy && errorarray(iters) < precisionGD
                                                    if printData
                                                        fprintf('GKV Parameters Converged to precision, Iteration No. %d.\n',iters);
                                                    end
                                                    endFlag = 1;
                                                end
                                            end

                                            if iters >= n_maxIterations
                                                if printData
                                                    fprintf('GKV Exceeded number of desired iterations, Loop %d.\n',j);
                                                end
                                                endFlag = 1;
                                            end

                                            if iters >= 1000 && (errorarray(iters) / initialError >= abandonThreshold)
                                                if printData
                                                    fprintf('GKV Error does not appear to converge, Loop %d. Continuing to next attempt.\n',j);
                                                end
                                                endFlag = 1;
                                            end

                                        end % End While

                                        if plotFit
                                            figure(fitPlot)
                                            scatter(x_fit_loop(:,1),y_fit_GD_loop,'bx')
                                            hold on
                                            plot(x_fit_loop(:,1),F_conv_loop(beta0_conv_loop,x_fit_loop),'b-')
                                            plot(x_fit_loop(:,1),F_conv_loop(newParams_voigt,x_fit_loop),'r-')
                                            title('Initial vs. Final Optimized Model [Voigt]')
                                            xlabel('Time [s]')
                                            ylabel('Action Integral Value')
                                            set(gca,'YScale','log')
                                            set(gca,'XScale','log')
                                            legend('Data','Initial','Final')
                                            hold off
                                        end

                                        if ~strcmp('',settleGD)
                                            switch settleGD
                                                case 'newton'
                                                    newParams_temp = newtonMethod(newParams_voigt,x_fit_loop,y_fit_GD_loop,scaling,100,'gkv',elasticSetting,fluidSetting,subref_loop,U_func_loop,dQdI,integralCalc,errorCalcVoigt,plotCost);
                                                    settleError = errorCalcVoigt(newParams_temp,x_fit_loop);

                                                case 'nelder-mead'
                                                    options = optimset('Display','none',...
                                                        'PlotFcns',[],...
                                                        'MaxFunEvals',n_maxIterations,...
                                                        'MaxIter',n_maxIterations,...
                                                        'TolFun',0,...
                                                        'TolX',0);

                                                    [newParams_temp,finalRes,eflag] = fminsearch(@(x) boundObjective(@(x) sum((F_conv_loop(x,x_fit_loop)-y_fit_GD_loop).^2),x,ub_loop,lb_loop), newParams_voigt, options);
                                                    settleError = errorCalcVoigt(newParams_temp,x_fit_loop);

                                                case 'gen-alg'
                                                    if license('test','GADS_Toolbox')
                                                        opts = optimoptions('ga',...
                                                            'FunctionTolerance',scaling,...
                                                            'UseVectorized',true,...
                                                            ...%'HybridFcn','fminunc',...
                                                            'PlotFcn',@gaplotbestf); % @gaplotgenealogy
                                                        [newParams_temp,finalRes,eflag,simoutput] = ga(@(x) sum((F_conv_loop(x,x_fit_loop)-y_fit_GD_loop).^2),length(newParams_voigt),[],[],[],[],lb_loop,ub_loop);
                                                    else
                                                        fprintf('You do not have the MATLAB Global Optimization Toolbox, and cannot use the Genetic Algorithm (ga) function.\nThe gradient descent results have NOT been fed forward.')
                                                    end
                                                    settleError = errorCalcVoigt(newParams_temp,x_fit_loop);

                                            end

                                            if settleError < errorCalcVoigt(newParams_voigt,x_fit_loop)
                                                if ( ~any((ub_loop - newParams_temp) < 0) ) && ( ~any((newParams_temp - lb_loop) < 0) )
                                                    newParams_voigt = newParams_temp;
                                                    if printData
                                                        fprintf('Settling Method (%s) improved the Fit for GKV\n',settleGD)
                                                    end
                                                else
                                                    if printData
                                                        fprintf('Settling Method (%s) gave negative parameters for GKV. Throwing away settled parameters.\n',settleGD)
                                                    end
                                                end
                                            end

                                            if plotFit
                                                figure(fitPlot)
                                                hold on
                                                plot(x_fit_loop(:,1),F_conv_loop(newParams_voigt,x_fit_loop),'g-')
                                                legend('Data','Initial','Final','Post-Settling')
                                                hold off
                                            end
                                        end

                                        beta_dist(:,j) = newParams_voigt;
                                        beta0_dist(:,j) = beta0_conv_loop;

                                        % Calculate Error Against Data Streams
                                        if forceResidual
                                            yModelIndentation = (1./alphaModel.*F_conv_loop(newParams_voigt,inputs_loop)).^(2/3);
                                            resnorm_dist(j) = sum(((yModelIndentation-ydata_voigt_loop).^2)./movvar(ydata_voigt_loop,3))./(length(ydata_voigt_loop)-length(newParams_voigt)); % Standard Error of Regression (S)
                                            voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',log_scale(yModelIndentation,x_fit_loop(:,1),dt,st)'];
                                        else
                                            resnorm_dist(j) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)
                                            voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',F_conv_wrapper_loop(newParams_voigt,x_fit_loop)'];
                                        end

                                        if timeOperation
                                            convTime(j) = Par.toc;
                                        end

                                        if plotCost
                                            figure(costPlot);
                                            clf;
                                        end

                                        % Create random starting point
                                        [beta0_maxwell_loop,tauInds,modulusInds] = makeRandomParams(beta0_maxwell,ub_randLimit_maxwell,lb_randLimit_maxwell,elasticSetting_maxwell,fluidSetting_maxwell,newInds_maxwell);

                                        weightArray = ones(size(beta0_maxwell_loop));
                                        if scaleParams
                                            numTau = length(tauInds);
                                            if strcmp(fluidSetting_maxwell,'y')
                                               numTau = numTau-1;
                                            end
                                            if i == 1
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_maxwell_loop_weight = logscaleParams(h_conv_maxwell_loop);
                                                        else
                                                            h_conv_maxwell_loop_weight = logscaleParams(Q_func_maxwell_loop);
                                                        end

                                                    otherwise
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(modulusInds) = 10.^((max(ub_randLimit_maxwell(modulusInds))-min(lb_randLimit_maxwell(modulusInds)))/2);
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);
                                                        else
                                                            h_conv_maxwell_loop_weight = weightParams(Q_func_maxwell_loop,weightArray);
                                                        end
                                                end
                                            else
                                                switch lower(paramScale)
                                                    case 'equal'
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_maxwell(intersect(find(~newInds_maxwell),modulusInds))))));
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);
                                                        else
                                                            h_conv_maxwell_loop_weight = weightParams(Q_func_maxwell_loop,weightArray);
                                                        end
                                                    case 'individual'
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(intersect(find(~newInds_maxwell),modulusInds)) = 10.^(floor(log10(beta0_maxwell(intersect(find(~newInds_maxwell),modulusInds)))));
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);
                                                        else
                                                            h_conv_maxwell_loop_weight = weightParams(Q_func_maxwell_loop,weightArray);
                                                        end
                                                    case 'logscale'
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_maxwell_loop_weight = logscaleParams(h_conv_maxwell_loop);
                                                        else
                                                            h_conv_maxwell_loop_weight = logscaleParams(Q_func_maxwell_loop);
                                                        end
                                                end
                                            end
                                        else
                                            if ~strcmp(GDmethod,'alrm-norm')
                                                h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);
                                            else
                                                h_conv_maxwell_loop_weight = weightParams(Q_func_maxwell_loop,weightArray);
                                            end
                                        end

                                        errorarray = NaN(1,n_maxIterations);
                                        iters = 0;
                                        endFlag = 0;
                                        perturbWait = 0;

                                        if plotCost
                                            legendEntries = {};
                                            for jj = 1:length(beta0_maxwell_loop)
                                                if ismember(jj,tauInds)
                                                    if strcmp(fluidSetting_maxwell,'y') && jj == length(beta0_maxwell_loop)
                                                        legendEntries = horzcat(legendEntries,sprintf('$\\phi_f$'));
                                                    else
                                                        legendEntries = horzcat(legendEntries,sprintf('$\\tau_%d$',find(tauInds==jj)));
                                                    end
                                                elseif ismember(jj,modulusInds)
                                                    if strcmp(elasticSetting_maxwell,'y') && jj == 1
                                                        legendEntries = horzcat(legendEntries,sprintf('$E_e$'));
                                                    else
                                                        legendEntries = horzcat(legendEntries,sprintf('$E_%d$',find(modulusInds==jj)));
                                                    end
                                                end
                                            end
                                        end

                                        % Gradient Descent
                                        if timeOperation
                                            Par.tic;
                                        end

                                        while ~endFlag

                                            iters = iters + 1;
                                            if iters > 1
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    currentParams = newParams_maxwell./weightArray;
                                                else
                                                    currentParams = log10(newParams_maxwell);
                                                end
                                            end
                                            currentDiff = zeros(size(beta0_maxwell_loop));
                                            if iters == 1
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    currentParams = beta0_maxwell_loop./weightArray;
                                                else
                                                    currentParams = log10(beta0_maxwell_loop);
                                                end
                                                alphaTemp = alpha;
                                                oldStep = zeros(size(currentDiff));
                                                gammaTemp = gamma;
                                                runningAvGradient = zeros(size(currentDiff));
                                                currentUpdate = zeros(size(currentDiff));
                                                runningAvUpdate = zeros(size(currentDiff));
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    initialError = errorCalcMaxwell(currentParams.*weightArray,maxwellInputs_loop); % Standard Error of Regression (S)
                                                else
                                                    initialError = errorCalcMaxwell(10.^currentParams,maxwellInputs_loop); % Standard Error of Regression (S)
                                                end

                                                if ~symbolicJacobian
                                                    dQdI = cell(size(beta0_maxwell_loop));
                                                    for jj = 1:length(beta0_maxwell_loop)
                                                        if jj == 1 && strcmp(elasticSetting_maxwell,'y')
                                                            dQdI{jj} = @(c,inputs) ones(size(inputs(:,1)));
        %                                                     dQdI{jj} = @(c,inputs) [1;zeros(length(inputs(:,1))-1,1)];
                                                        elseif any(jj == modulusInds)
                                                            tempFunc = sprintf('@(c,inputs) exp( (-inputs(:,1))./c(%d) )',jj+1);
                                                            dQdI{jj} = str2func(tempFunc);
                                                        elseif any(jj == tauInds)
                                                            tempFunc = sprintf('@(c,inputs) (c(%d).*inputs(:,1)./(c(%d).^2)).*exp( (-inputs(:,1))./c(%d) )',jj-1,jj,jj);
                                                            dQdI{jj} = str2func(tempFunc);
                                                        else
                                                            error('Error creating the derivative for term %d!',jj)
                                                        end
                                                    end
                                                else
                                                    [J_sym,H_sym,dQdI] = getObjectiveDerivatives(beta0_maxwell_loop,'gm',elasticSetting_maxwell,fluidSetting_maxwell);
                                                end
                                            end

                                            % Calculate gradients inside integral term
                                            for jj = 1:length(beta0_maxwell_loop)
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        if ~normalizeGD
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = (1/length(y_fit_maxwell_GD_loop)).*sum(integralCalc(10.^currentParams,inputsForce,dQdI{jj})-y_fit_maxwell_GD_loop); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = (2/length(y_fit_maxwell_GD_loop)).*sum((h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*integralCalc(10.^currentParams,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = (2/length(y_fit_maxwell_GD_loop)).*sum((h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*jacobianCalc(10.^currentParams,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        else
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = sum((1+100.*inputsForce(:,1).^2).*(integralCalc(10.^currentParams,inputsForce,dQdI{jj})-y_fit_maxwell_GD_loop)); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce(:,1).^2).*(h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*integralCalc(10.^currentParams,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce(:,1).^2).*(h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*jacobianCalc(10.^currentParams,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        end

                                                    otherwise
                                                        if ~normalizeGD
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = (1/length(y_fit_maxwell_GD_loop)).*sum(integralCalc(currentParams.*weightArray,inputsForce,dQdI{jj})-y_fit_maxwell_GD_loop); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = (2/length(y_fit_maxwell_GD_loop)).*sum((h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*integralCalc(currentParams.*weightArray,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = (2/length(y_fit_maxwell_GD_loop)).*sum((h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*jacobianCalc(currentParams.*weightArray,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        else
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = sum((1+100.*inputsForce(:,1).^2).*(integralCalc(currentParams.*weightArray,inputsForce,dQdI{jj})-y_fit_maxwell_GD_loop)); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce(:,1).^2).*(h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*integralCalc(currentParams.*weightArray,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce(:,1).^2).*(h_conv_maxwell_loop_weight(currentParams,inputsForce)-y_fit_maxwell_GD_loop).*jacobianCalc(currentParams.*weightArray,inputsForce,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        end
                                                end
                                            end

                                            % Update current parameters
                                            switch GDmethod
                                                case 'standard'
                                                    newParams_maxwell = currentParams - alphaTemp.*currentDiff;
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_maxwell = 10.^newParams_maxwell;
                                                        otherwise
                                                            newParams_maxwell = newParams_maxwell.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            newParams_maxwell = currentParams + rand(size(currentParams)).*currentParams;
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_maxwell) < 0) ) || ( any((newParams_maxwell - lb_loop_maxwell) < 0) )
                                                        newParams_maxwell((ub_loop_maxwell - newParams_maxwell) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_maxwell) < 0);
                                                        newParams_maxwell((newParams_maxwell - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_maxwell - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)
                                                    alphaarray(:,iters) = -alphaTemp.*currentDiff;

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,alpha,'rx')
                                                        title('Parameter Update vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Parameter Update')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case 'sgd'
                                                    alphaarray(:,iters) = -alphaTemp.*(currentDiff.*(sqrt(3)*(2*rand(size(currentDiff))-1)));
                                                    newParams_maxwell = currentParams + alphaarray(:,iters)';
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_maxwell = 10.^newParams_maxwell;
                                                        otherwise
                                                            newParams_maxwell = newParams_maxwell.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            newParams_maxwell = currentParams + rand(size(currentParams)).*currentParams;
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_maxwell) < 0) ) || ( any((newParams_maxwell - lb_loop_maxwell) < 0) )
                                                        newParams_maxwell((ub_loop_maxwell - newParams_maxwell) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_maxwell) < 0);
                                                        newParams_maxwell((newParams_maxwell - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_maxwell - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Parameter Update vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Parameter Update')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case {'alrm','alrm-norm'}
                                                    newParams_maxwell = currentParams + (gammaTemp.*oldStep - (1-gammaTemp).*alphaTemp.*currentDiff);
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_maxwell = 10.^newParams_maxwell;
                                                        otherwise
                                                            newParams_maxwell = newParams_maxwell.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            switch perturbMethod
                                                                case 'params'
                                                                    newParams_maxwell = currentParams + rand(size(currentParams)).*currentParams;
                                                                case 'LR'
                                                                    alphaTemp = alpha;
                                                            end
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_maxwell) < 0) ) || ( any((newParams_maxwell - lb_loop_maxwell) < 0) )
                                                        newParams_maxwell((ub_loop_maxwell - newParams_maxwell) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_maxwell) < 0);
                                                        newParams_maxwell((newParams_maxwell - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_maxwell - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)

                                                    % Test if a valid step was taken
                                                    if iters == 1
                                                        errorChange = errorarray(iters) / initialError; % Error Ratio
                                                    else
                                                        errorChange = errorarray(iters) / errorarray(iters-1); % Error Ratio
                                                    end

                                                    if errorChange < zeta
                                                        alphaTemp = alphaTemp * eta; % Increase learning rate
                                                        gammaTemp = gamma;
                                                        oldStep = (gammaTemp.*oldStep - (1-gammaTemp).*alphaTemp.*currentDiff);
                                                    else
                                                        alphaTemp = alphaTemp * rho; % Decrease learning rate
                                                        newParams_maxwell = currentParams; % Throw away parameter updates
                                                        switch lower(paramScale)
                                                            case 'logscale'
                                                                newParams_maxwell = 10.^newParams_maxwell;
                                                            otherwise
                                                                newParams_maxwell = newParams_maxwell.*weightArray;
                                                        end
                                                        errorarray(iters) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)
                                                        gammaTemp = 0;
                                                    end
                                                    alphaarray(iters) = alphaTemp;

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        scatter((1:iters),alphaarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,alpha,'rx')
                                                        title('Learning Parameter vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Learning Parameter')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case 'adadelta'
                                                    runningAvGradient = rho_ada.*(runningAvGradient)+(1-rho_ada).*(currentDiff.^2);
                                                    currentUpdate = -(sqrt(runningAvUpdate+epsilon_adadelta)./sqrt(runningAvGradient+epsilon_adadelta)).*currentDiff;
                                                    runningAvUpdate = rho_ada.*(runningAvUpdate)+(1-rho_ada).*(currentUpdate.^2);
                                                    newParams_maxwell = currentParams + currentUpdate;
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_maxwell = 10.^newParams_maxwell;
                                                        otherwise
                                                            newParams_maxwell = newParams_maxwell.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            switch perturbMethod
                                                                case 'params'
                                                                    newParams_maxwell = currentParams + rand(size(currentParams)).*currentParams;
                                                                case 'LR'
                                                                    runningAvUpdate = (runningAvUpdate./abs(runningAvUpdate)).*abs(runningAvUpdate).^(rand(size(runningAvUpdate))+1);
                                                            end
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_maxwell) < 0) ) || ( any((newParams_maxwell - lb_loop_maxwell) < 0) )
                                                        newParams_maxwell((ub_loop_maxwell - newParams_maxwell) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_maxwell) < 0);
                                                        newParams_maxwell((newParams_maxwell - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_maxwell - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)
                                                    alphaarray(:,iters) = (sqrt(runningAvUpdate)./sqrt(runningAvGradient+epsilon_adadelta));

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Learning Parameter vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Learning Parameter')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                            end

                                            perturbWait = perturbWait + 1;

                                            % Check for stop conditions
                                            if iters > 1
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        temp = 10.^currentParams;
                                                    otherwise
                                                        temp = currentParams.*weightArray;
                                                end
                                                if max(abs((temp-newParams_maxwell)./temp)) < accuracy && errorarray(iters) < precisionGD
                                                    if printData
                                                        fprintf('GM Parameters Converged to precision, Iteration No. %d.\n',iters);
                                                    end
                                                    endFlag = 1;
                                                end
                                            end

                                            if iters >= n_maxIterations
                                                if printData
                                                    fprintf('GM Exceeded number of desired iterations, Loop %d.\n',j);
                                                end
                                                endFlag = 1;
                                            end

                                            if iters >= 1000 && (errorarray(iters) / initialError >= abandonThreshold)
                                                if printData
                                                    fprintf('GM Error does not appear to converge, Loop %d. Continuing to next attempt.\n',j);
                                                end
                                                endFlag = 1;
                                            end

                                        end % End While

                                        if plotFit
                                            figure(fitPlot)
                                            scatter(maxwellInputs_loop(:,1),y_fit_maxwell_GD_loop,'bx')
                                            hold on
                                            plot(maxwellInputs_loop(:,1),h_conv_maxwell_loop(beta0_maxwell_loop,maxwellInputs_loop),'b-')
                                            plot(maxwellInputs_loop(:,1),h_conv_maxwell_loop(newParams_maxwell,maxwellInputs_loop),'r-')
                                            title('Initial vs. Final Optimized Model [Maxwell]')
                                            xlabel('Time [s]')
                                            ylabel('Action Integral Value')
                                            set(gca,'YScale','log')
                                            set(gca,'XScale','log')
                                            legend('Data','Initial','Final')
                                            hold off
                                        end

                                        if ~strcmp('',settleGD)
                                            switch settleGD
                                                case 'newton'
                                                    newParams_temp = newtonMethod(newParams_maxwell,maxwellInputs_loop,y_fit_maxwell_GD_loop,scaling,100,'gm',elasticSetting_maxwell,fluidSetting_maxwell,subref_loop,Q_func_maxwell_loop,dQdI,integralCalc,errorCalcMaxwell,plotCost);
                                                    settleError = errorCalcMaxwell(newParams_temp,maxwellInputs_loop);

                                                case 'nelder-mead'
                                                    options = optimset('Display','none',...
                                                        'PlotFcns',[],...
                                                        'MaxFunEvals',n_maxIterations,...
                                                        'MaxIter',n_maxIterations,...
                                                        'TolFun',0,...
                                                        'TolX',0);

                                                    [newParams_temp,finalRes,eflag] = fminsearch(@(x) boundObjective(@(x) sum((h_conv_maxwell_loop(x,maxwellInputs_loop)-y_fit_maxwell_GD_loop).^2),x,ub_loop_maxwell,lb_loop_maxwell), newParams_maxwell, options);
                                                    settleError = errorCalcMaxwell(newParams_temp,maxwellInputs_loop);

                                                case 'gen-alg'
                                                    if license('test','GADS_Toolbox')
                                                        opts = optimoptions('ga',...
                                                            'FunctionTolerance',scaling,...
                                                            'UseVectorized',true,...
                                                            ...%'HybridFcn','fminunc',...
                                                            'PlotFcn',@gaplotbestf); % @gaplotgenealogy
                                                        [newParams_temp,finalRes,eflag,simoutput] = ga(@(x) sum((h_conv_maxwell_loop(x,maxwellInputs_loop)-y_fit_maxwell_GD_loop).^2),length(newParams_maxwell),[],[],[],[],lb_loop_maxwell,ub_loop_maxwell);
                                                    else
                                                        fprintf('You do not have the MATLAB Global Optimization Toolbox, and cannot use the Genetic Algorithm (ga) function.\nThe gradient descent results have NOT been fed forward.')
                                                    end
                                                    settleError = errorCalcMaxwell(newParams_temp,maxwellInputs_loop);

                                            end

                                            if settleError < errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop)
                                                if ( ~any((ub_loop_maxwell - newParams_temp) < 0) ) && ( ~any((newParams_temp - lb_loop_maxwell) < 0) )
                                                    newParams_maxwell = newParams_temp;
                                                    if printData
                                                        fprintf('Settling Method (%s) improved the Fit for GM\n',settleGD)
                                                    end
                                                else
                                                    if printData
                                                        fprintf('Settling Method (%s) gave negative parameters for GKV. Throwing away settled parameters.\n',settleGD)
                                                    end
                                                end
                                            end

                                            if plotFit
                                                figure(fitPlot)
                                                hold on
                                                plot(maxwellInputs_loop(:,1),h_conv_maxwell_loop(newParams_maxwell,maxwellInputs_loop),'g-')
                                                legend('Data','Initial','Final','Post-Settling')
                                                hold off
                                            end
                                        end

                                        beta_dist_maxwell(:,j) = newParams_maxwell;
                                        beta0_dist_maxwell(:,j) = beta0_maxwell_loop;

                                        % Calculate Error Against Data Streams
                                        if forceResidual
                                            if strcmp(elasticSetting_maxwell,'y')
                                                yModelForce = alphaModel.*(newParams_maxwell(1).*(inputsForce(:,2))+subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                            else
                                                yModelForce = alphaModel.*(subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                            end
                                            resnorm_dist_maxwell(j) = sum(((yModelForce-ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(newParams_maxwell)); % Standard Error of Regression (S)
                                            maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',log_scale(yModelForce,maxwellInputs_loop(:,1),dt,st)'];
                                        else
                                            resnorm_dist_maxwell(j) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)
                                            maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',h_conv_maxwell_wrapper_loop(newParams_maxwell,maxwellInputs_loop)'];
                                        end

                                        if timeOperation
                                            maxwellTime(j) = Par.toc;
                                        end

                                        if plotCost
                                            figure(costPlot);
                                            clf;
                                        end

                                        % Create random starting point
                                        [beta0_schapery_loop,tauInds,modulusInds] = makeRandomParams(beta0_schapery,ub_randLimit_schapery,lb_randLimit_schapery,elasticSetting_maxwell,fluidSetting_maxwell,newInds_schapery);

                                        weightArray = ones(size(beta0_schapery_loop));
                                        if scaleParams
                                            numTau = length(tauInds);
                                            if strcmp(fluidSetting_maxwell,'y')
                                               numTau = numTau-1;
                                            end
                                            if i == 1
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_schapery_area_loop_weight = logscaleParams(h_conv_schapery_area_loop);
                                                        else
                                                            h_conv_schapery_area_loop_weight = logscaleParams(Q_func_loop);
                                                        end

                                                    otherwise
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(modulusInds) = 10.^((max(ub_randLimit_schapery(modulusInds))-min(lb_randLimit_schapery(modulusInds)))/2);
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);
                                                        else
                                                            h_conv_schapery_area_loop_weight = weightParams(Q_func_loop,weightArray);
                                                        end
                                                end
                                            else
                                                switch lower(paramScale)
                                                    case 'equal'
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_schapery(intersect(find(~newInds_schapery),modulusInds))))));
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);
                                                        else
                                                            h_conv_schapery_area_loop_weight = weightParams(Q_func_loop,weightArray);
                                                        end
                                                    case 'individual'
                                                        weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                        weightArray(intersect(find(~newInds_schapery),modulusInds)) = 10.^(floor(log10(beta0_schapery(intersect(find(~newInds_schapery),modulusInds)))));
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);
                                                        else
                                                            h_conv_schapery_area_loop_weight = weightParams(Q_func_loop,weightArray);
                                                        end
                                                    case 'logscale'
                                                        if ~strcmp(GDmethod,'alrm-norm')
                                                            h_conv_schapery_area_loop_weight = logscaleParams(h_conv_schapery_area_loop);
                                                        else
                                                            h_conv_schapery_area_loop_weight = logscaleParams(Q_func_loop);
                                                        end
                                                end
                                            end
                                        else
                                            if ~strcmp(GDmethod,'alrm-norm')
                                                h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);
                                            else
                                                h_conv_schapery_area_loop_weight = weightParams(Q_func_loop,weightArray);
                                            end
                                        end

                                        errorarray = NaN(1,n_maxIterations);
                                        iters = 0;
                                        endFlag = 0;
                                        perturbWait = 0;

                                        if plotCost
                                            legendEntries = {};
                                            for jj = 1:length(beta0_maxwell_loop)
                                                if ismember(jj,tauInds)
                                                    if strcmp(fluidSetting_maxwell,'y') && jj == length(beta0_maxwell_loop)
                                                        legendEntries = horzcat(legendEntries,sprintf('$\\phi_f$'));
                                                    else
                                                        legendEntries = horzcat(legendEntries,sprintf('$\\tau_%d$',find(tauInds==jj)));
                                                    end
                                                elseif ismember(jj,modulusInds)
                                                    if strcmp(elasticSetting_maxwell,'y') && jj == 1
                                                        legendEntries = horzcat(legendEntries,sprintf('$E_e$'));
                                                    else
                                                        legendEntries = horzcat(legendEntries,sprintf('$E_%d$',find(modulusInds==jj)));
                                                    end
                                                end
                                            end
                                        end

                                        % Gradient Descent
                                        if timeOperation
                                            Par.tic;
                                        end

                                        while ~endFlag

                                            iters = iters + 1;
                                            if iters > 1
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    currentParams = newParams_schapery./weightArray;
                                                else
                                                    currentParams = log10(newParams_schapery);
                                                end
                                            end
                                            currentDiff = zeros(size(beta0_schapery_loop));
                                            if iters == 1
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    currentParams = beta0_schapery_loop./weightArray;
                                                else
                                                    currentParams = log10(beta0_schapery_loop);
                                                end
                                                alphaTemp = alphaS;
                                                oldStep = zeros(size(currentDiff));
                                                gammaTemp = gammaS;
                                                runningAvGradient = zeros(size(currentDiff));
                                                currentUpdate = zeros(size(currentDiff));
                                                runningAvUpdate = zeros(size(currentDiff));
                                                if ~strcmp(lower(paramScale),'logscale')
                                                    initialError = errorCalcSchapery(currentParams.*weightArray,x_fit_area_loop); % Standard Error of Regression (S)
                                                else
                                                    initialError = errorCalcSchapery(10.^currentParams,x_fit_area_loop); % Standard Error of Regression (S)
                                                end
                                                if ~symbolicJacobian
                                                    dQdI = cell(size(beta0_schapery_loop));
                                                    for jj = 1:length(beta0_schapery_loop)
                                                        if jj == 1 && strcmp(elasticSetting_maxwell,'y')
                                                            dQdI{jj} = @(c,inputs) ones(size(inputs(:,1)));
        %                                                     dQdI{jj} = @(c,inputs) [1;zeros(length(inputs(:,1))-1,1)];
                                                        elseif any(jj == modulusInds)
                                                            tempFunc = sprintf('@(c,inputs) exp( (-inputs(:,1))./c(%d) )',jj+1);
                                                            dQdI{jj} = str2func(tempFunc);
                                                        elseif any(jj == tauInds)
                                                            tempFunc = sprintf('@(c,inputs) (c(%d).*inputs(:,1)./(c(%d).^2)).*exp( (-inputs(:,1))./c(%d) )',jj-1,jj,jj);
                                                            dQdI{jj} = str2func(tempFunc);
                                                        else
                                                            error('Error creating the derivative for term %d!',jj)
                                                        end
                                                    end
                                                else
                                                    [J_sym,H_sym,dQdI] = getObjectiveDerivatives(beta0_schapery_loop,'gm',elasticSetting_maxwell,fluidSetting_maxwell);
                                                end
                                            end

                                            % Calculate gradients inside integral term
                                            for jj = 1:length(beta0_schapery_loop)
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        if ~normalizeGD
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = (1/length(ydata_loop)).*sum(integralCalc(10.^currentParams,inputsForce_schapery,dQdI{jj})-ydata_loop); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = (2/length(ydata_loop)).*sum((h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*integralCalc(10.^currentParams,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = (2/length(ydata_loop)).*sum((h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*jacobianCalc(10.^currentParams,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        else
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = sum((1+100.*inputsForce_schapery(:,1).^2).*(integralCalc(10.^currentParams,inputsForce_schapery,dQdI{jj})-ydata_loop)); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce_schapery(:,1).^2).*(h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*integralCalc(10.^currentParams,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce_schapery(:,1).^2).*(h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*jacobianCalc(10.^currentParams,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        end

                                                    otherwise
                                                        if ~normalizeGD
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = (1/length(ydata_loop)).*sum(integralCalc(currentParams.*weightArray,inputsForce_schapery,dQdI{jj})-ydata_loop); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = (2/length(ydata_loop)).*sum((h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*integralCalc(currentParams.*weightArray,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = (2/length(ydata_loop)).*sum((h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*jacobianCalc(currentParams.*weightArray,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        else
                                                            switch GDmethod
                                                                case 'standard'
                                                                    currentDiff(jj) = sum((1+100.*inputsForce_schapery(:,1).^2).*(integralCalc(currentParams.*weightArray,inputsForce_schapery,dQdI{jj})-ydata_loop)); % Derivative of Residual
                                                                case {'alrm','adadelta','sgd'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce_schapery(:,1).^2).*(h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*integralCalc(currentParams.*weightArray,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                                case{'alrm-norm'}
                                                                    currentDiff(jj) = 2*sum((1+100.*inputsForce_schapery(:,1).^2).*(h_conv_schapery_area_loop_weight(currentParams,inputsForce_schapery)-ydata_loop).*jacobianCalc(currentParams.*weightArray,inputsForce_schapery,dQdI{jj})); % Hackerearth Article for Batch GD
                                                            end
                                                        end
                                                end
                                            end

                                            % Update current parameters
                                            switch GDmethod
                                                case 'standard'
                                                    newParams_schapery = currentParams - alphaTemp.*currentDiff;
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_schapery = 10.^newParams_schapery;
                                                        otherwise
                                                            newParams_schapery = newParams_schapery.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            newParams_schapery = currentParams + rand(size(currentParams)).*currentParams;
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_schapery) < 0) ) || ( any((newParams_schapery - lb_loop_maxwell) < 0) )
                                                        newParams_schapery((ub_loop_maxwell - newParams_schapery) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_schapery) < 0);
                                                        newParams_schapery((newParams_schapery - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_schapery - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)
                                                    alphaarray(:,iters) = -alphaTemp.*currentDiff;

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Schapery-Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Parameter Update vs. Iteration No. [Schapery-Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Parameter Update')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case 'sgd'
                                                    alphaarray(:,iters) = -alphaTemp.*(currentDiff.*(sqrt(3)*(2*rand(size(currentDiff))-1)));
                                                    newParams_schapery = currentParams + alphaarray(:,iters)';
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_schapery = 10.^newParams_schapery;
                                                        otherwise
                                                            newParams_schapery = newParams_schapery.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            newParams_schapery = currentParams + rand(size(currentParams)).*currentParams;
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_schapery) < 0) ) || ( any((newParams_schapery - lb_loop_maxwell) < 0) )
                                                        newParams_schapery((ub_loop_maxwell - newParams_schapery) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_schapery) < 0);
                                                        newParams_schapery((newParams_schapery - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_schapery - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Parameter Update vs. Iteration No. [Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Parameter Update')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case {'alrm','alrm-norm'}
                                                    newParams_schapery = currentParams + (gammaTemp.*oldStep - (1-gammaTemp).*alphaTemp.*currentDiff);
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_schapery = 10.^newParams_schapery;
                                                        otherwise
                                                            newParams_schapery = newParams_schapery.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            switch perturbMethod
                                                                case 'params'
                                                                    newParams_schapery = currentParams + rand(size(currentParams)).*currentParams;
                                                                case 'LR'
                                                                    alphaTemp = alphaS;
                                                            end
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_schapery) < 0) ) || ( any((newParams_schapery - lb_loop_maxwell) < 0) )
                                                        newParams_schapery((ub_loop_maxwell - newParams_schapery) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_schapery) < 0);
                                                        newParams_schapery((newParams_schapery - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_schapery - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)

                                                    % Test if a valid step was taken
                                                    if iters == 1
                                                        errorChange = errorarray(iters) / initialError; % Error Ratio
                                                    else
                                                        errorChange = errorarray(iters) / errorarray(iters-1); % Error Ratio
                                                    end

                                                    if errorChange < zeta
                                                        alphaTemp = alphaTemp * eta; % Increase learning rate
                                                        gammaTemp = gammaS;
                                                        oldStep = (gammaTemp.*oldStep - (1-gammaTemp).*alphaTemp.*currentDiff);
                                                    else
                                                        alphaTemp = alphaTemp * rho; % Decrease learning rate
                                                        newParams_schapery = currentParams; % Throw away parameter updates
                                                        switch lower(paramScale)
                                                            case 'logscale'
                                                                newParams_schapery = 10.^newParams_schapery;
                                                            otherwise
                                                                newParams_schapery = newParams_schapery.*weightArray;
                                                        end
                                                        errorarray(iters) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)
                                                        gammaTemp = 0;
                                                    end
                                                    alphaarray(iters) = alphaTemp;

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Schapery-Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        scatter((1:iters),alphaarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,alphaS,'rx')
                                                        title('Learning Parameter vs. Iteration No. [Schapery-Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Learning Parameter')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                    end

                                                case 'adadelta'
                                                    runningAvGradient = rho_adaS.*(runningAvGradient)+(1-rho_adaS).*(currentDiff.^2);
                                                    currentUpdate = -(sqrt(runningAvUpdate+epsilon_adadelta)./sqrt(runningAvGradient+epsilon_adadelta)).*currentDiff;
                                                    runningAvUpdate = rho_adaS.*(runningAvUpdate)+(1-rho_adaS).*(currentUpdate.^2);
                                                    newParams_schapery = currentParams + currentUpdate;
                                                    switch lower(paramScale)
                                                        case 'logscale'
                                                            newParams_schapery = 10.^newParams_schapery;
                                                        otherwise
                                                            newParams_schapery = newParams_schapery.*weightArray;
                                                    end

                                                    % Check for perturbation condition
                                                    if iters > perturbWindow && perturbGD && (perturbWait > perturbWindow)
                                                        if max(abs(errorarray(iters-perturbWindow)-mean(errorarray((iters-perturbWindow):iters-1)))) < perturbBand*mean(errorarray((iters-perturbWindow):(iters-1)))
                                                            switch perturbMethod
                                                                case 'params'
                                                                    newParams_schapery = currentParams + rand(size(currentParams)).*currentParams;
                                                                case 'LR'
                                                                    runningAvUpdate = (runningAvUpdate./abs(runningAvUpdate)).*abs(runningAvUpdate).^(rand(size(runningAvUpdate))+1);
                                                            end
                                                            perturbWait = 0;
                                                        end
                                                    end

                                                    if ( any((ub_loop_maxwell - newParams_schapery) < 0) ) || ( any((newParams_schapery - lb_loop_maxwell) < 0) )
                                                        newParams_schapery((ub_loop_maxwell - newParams_schapery) < 0) = ub_loop_maxwell((ub_loop_maxwell - newParams_schapery) < 0);
                                                        newParams_schapery((newParams_schapery - lb_loop_maxwell) < 0) = lb_loop_maxwell((newParams_schapery - lb_loop_maxwell) < 0);
                                                    end
                                                    errorarray(iters) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)
                                                    alphaarray(:,iters) = (sqrt(runningAvUpdate)./sqrt(runningAvGradient+epsilon_adadelta));

                                                    if plotCost && ~mod(iters,n_plot)
                                                        figure(costPlot)
                                                        subplot(2,1,1)
                                                        scatter((1:iters),errorarray(1:iters),'bo')
                                                        hold on
                                                        scatter(0,initialError,'rx')
                                                        title('Cost vs. Iteration No. [Schapery-Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Cost')
                                                        set(gca,'YScale','log')
                                                        hold off
                                                        subplot(2,1,2)
                                                        plot((1:iters),alphaarray(:,1:iters))
                                                        hold on
                                                        scatter(0,0,'rx')
                                                        title('Learning Parameter vs. Iteration No. [Schapery-Maxwell]')
                                                        xlabel('Iteration No.')
                                                        ylabel('Learning Parameter')
                                                        legend(legendEntries,'interpreter','latex')
            %                                             set(gca,'YScale','log')
                                                        hold off
                                                    end

                                            end

                                            perturbWait = perturbWait + 1;

                                            if iters > 1
                                                switch lower(paramScale)
                                                    case 'logscale'
                                                        temp = 10.^currentParams;
                                                    otherwise
                                                        temp = currentParams.*weightArray;
                                                end
                                                if max(abs((temp-newParams_schapery)./temp)) < accuracy && errorarray(iters) < precisionGD
                                                    if printData
                                                        fprintf('SGM Parameters Converged to precision, Iteration No. %d.\n',iters);
                                                    end
                                                    endFlag = 1;
                                                end
                                            end

                                            if iters >= n_maxIterations
                                                if printData
                                                    fprintf('SGM Exceeded number of desired iterations, Loop %d.\n',j);
                                                end
                                                endFlag = 1;
                                            end

                                            if iters >= 1000 && (errorarray(iters) / initialError >= abandonThreshold)
                                                if printData
                                                    fprintf('SGM Error does not appear to converge, Loop %d. Continuing to next attempt.\n',j);
                                                end
                                                endFlag = 1;
                                            end

                                        end % End While

                                        if plotFit
                                            figure(fitPlot)
                                            scatter(x_fit_area_loop(:,1),y_fit_schapery_GD_loop,'bx')
                                            hold on
                                            plot(x_fit_area_loop(:,1),h_conv_schapery_area_loop(beta0_schapery_loop,x_fit_area_loop),'b-')
                                            plot(x_fit_area_loop(:,1),h_conv_schapery_area_loop(newParams_schapery,x_fit_area_loop),'r-')
                                            title('Initial vs. Final Optimized Model [Schapery-Maxwell]')
                                            xlabel('Time [s]')
                                            ylabel('Force [N]')
                                            set(gca,'YScale','log')
                                            set(gca,'XScale','log')
                                            legend('Data','Initial','Final')
                                            hold off
                                        end

                                        if ~strcmp('',settleGD)
                                            switch settleGD
                                                case 'newton'
                                                    newParams_temp = newtonMethod(newParams_schapery,x_fit_area_loop,y_fit_schapery_GD_loop,scaling,100,'gm',elasticSetting_maxwell,fluidSetting_maxwell,subref_loop,U_func_loop,dQdI,integralCalc,errorCalcSchapery,plotCost);
                                                    settleError = errorCalcSchapery(newParams_temp,x_fit_area_loop);

                                                case 'nelder-mead'
                                                    options = optimset('Display','none',...
                                                        'PlotFcns',[],...
                                                        'MaxFunEvals',n_maxIterations,...
                                                        'MaxIter',n_maxIterations,...
                                                        'TolFun',0,...
                                                        'TolX',0);

                                                    [newParams_temp,finalRes,eflag] = fminsearch(@(x) boundObjective(@(x) sum((h_conv_schapery_area_loop(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2),x,ub_loop_maxwell,lb_loop_maxwell), newParams_schapery, options);
                                                    settleError = errorCalcSchapery(newParams_temp,x_fit_area_loop);

                                                case 'gen-alg'
                                                    if license('test','GADS_Toolbox')
                                                        opts = optimoptions('ga',...
                                                            'FunctionTolerance',scaling,...
                                                            'UseVectorized',true,...
                                                            ...%'HybridFcn','fminunc',...
                                                            'PlotFcn',@gaplotbestf); % @gaplotgenealogy
                                                        [newParams_temp,finalRes,eflag,simoutput] = ga(@(x) sum((h_conv_schapery_area_loop(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2),length(newParams_schapery),[],[],[],[],lb_loop_maxwell,ub_loop_maxwell);
                                                    else
                                                        fprintf('You do not have the MATLAB Global Optimization Toolbox, and cannot use the Genetic Algorithm (ga) function.\nThe gradient descent results have NOT been fed forward.')
                                                    end
                                                    settleError = errorCalcSchapery(newParams_temp,x_fit_area_loop);

                                            end

                                            if settleError < errorCalcSchapery(newParams_schapery,x_fit_area_loop)
                                                if ( ~any((ub_loop_maxwell - newParams_temp) < 0) ) && ( ~any((newParams_temp - lb_loop_maxwell) < 0) )
                                                    newParams_schapery = newParams_temp;
                                                    if printData
                                                        fprintf('Settling Method (%s) improved the Fit for SGM\n',settleGD)
                                                    end
                                                else
                                                    if printData
                                                        fprintf('Settling Method (%s) gave negative parameters for GKV. Throwing away settled parameters.\n',settleGD)
                                                    end
                                                end
                                            end

                                            if plotFit
                                                figure(fitPlot)
                                                hold on
                                                plot(x_fit_area_loop(:,1),h_conv_schapery_area_loop(newParams_schapery,x_fit_area_loop),'g-')
                                                legend('Data','Initial','Final','Post-Settling')
                                                hold off
                                            end
                                        end

                                        beta_dist_schapery_area(:,j) = newParams_schapery;
                                        beta0_dist_schapery_area(:,j) = beta0_schapery_loop;

                                        % Calculate Error Against Data Streams
                                        if forceResidual
                                            yModelForce = h_conv_schapery_area_loop(newParams_schapery,inputsForce_schapery);
                                            resnorm_dist_schapery_area(j) = sum(((yModelForce-ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(newParams_schapery)); % Standard Error of Regression (S)
                                            schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',log_scale(yModelForce,x_fit_area_loop(:,1),dt,st)'];
                                        else
                                            resnorm_dist_schapery_area(j) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)
                                            schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',h_conv_schapery_area_wrapper_loop(newParams_schapery,x_fit_area_loop)'];
                                        end

                                        if timeOperation
                                            schaperyTime(j) = Par.toc;
                                        end

                                        if timeOperation
                                            % Track memory
                                            %ramused(j) = userMem.MemUsedMATLAB;
                                            ramused(j) = 0;
                                            %ramav(j) = systemMem.PhysicalMemory.Available;
                                            ramav(j) = 0;
                                            ramtime(j) = now;
                                        end

                                        % Placeholder for the linear model since we
                                        % don't care about it.
                                        beta_dist_linear(:,j) = ones(size(beta0_conv_loop));
                                        beta0_dist_linear(:,j) = ones(size(beta0_conv_loop));
                                        resnorm_dist_linear(j) = 1;

                                    end % End Parfor

                                    % Check the multiGrid status
                                    if multiGridGD
                                        if multiGrid_iters < multiGrid_hardCap
                                            bestRes = [min(resnorm_dist) ...
                                                min(resnorm_dist_maxwell)];
                                            if any(bestRes < precisionMultigrid)
                                                fprintf('\nMultiGrid resulted in convergence after %d grid(s)!', multiGrid_iters);
                                                multiGridFlag = 0;
                                            else
                                                % Tighten random search bounds by
                                                % creating error band around the
                                                % best parameter sets.
                                                resnorm_dist_temp = resnorm_dist;
                                                resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                                                [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                                                bestParams = beta_dist(:,idx);
                                                stdBand = std(beta_dist,0,2);
                                                ub_randLimit_multigrid = (bestParams + stdBand);
                                                lb_randLimit_multigrid = (bestParams - stdBand);

                                                if any(ub_randLimit_multigrid > ub_loop','all') || any(lb_randLimit_multigrid < lb_loop','all') || any(isnan(bestParams),'all') || any(isinf(bestParams),'all')
                                                   badInds = logical(sum([ub_randLimit_multigrid > ub_loop',isnan(bestParams),isinf(bestParams)],2));
                                                   ub_randLimit_multigrid(badInds) = ub_loop(badInds)';
                                                   badInds = logical(sum([lb_randLimit_multigrid < lb_loop',isnan(bestParams),isinf(bestParams)],2));
                                                   lb_randLimit_multigrid(badInds) = lb_loop(badInds)';
                                                end

                                                if strcmp(elasticSetting,'y')
                                                    if strcmp(fluidSetting,'y')
                                                        modulusInds = horzcat(1,(2:2:length(bestParams)-1));
                                                    else
                                                        modulusInds = horzcat(1,(2:2:length(bestParams)));
                                                    end
                                                else
                                                    if strcmp(fluidSetting,'y')
                                                        modulusInds = (1:2:length(bestParams)-1);
                                                    else
                                                        modulusInds = (1:2:length(bestParams));
                                                    end
                                                end

                                                ub_randLimit_multigrid(modulusInds) = log10(ub_randLimit_multigrid(modulusInds));
                                                lb_randLimit_multigrid(modulusInds) = log10(lb_randLimit_multigrid(modulusInds));

                                                resnorm_dist_temp = resnorm_dist_maxwell;
                                                resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                                                [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                                                bestParams = beta_dist_maxwell(:,idx);
                                                stdBand = std(beta_dist_maxwell,0,2);
                                                ub_randLimit_maxwell_multigrid = (bestParams + stdBand);
                                                lb_randLimit_maxwell_multigrid = (bestParams - stdBand);

                                                if any(ub_randLimit_maxwell_multigrid > ub_loop_maxwell','all') || any(lb_randLimit_maxwell_multigrid < lb_loop_maxwell','all') || any(isnan(bestParams),'all') || any(isinf(bestParams),'all')
                                                   badInds = logical(sum([ub_randLimit_maxwell_multigrid > ub_loop_maxwell',isnan(bestParams),isinf(bestParams)],2));
                                                   ub_randLimit_maxwell_multigrid(badInds) = ub_loop_maxwell(badInds)';
                                                   badInds = logical(sum([lb_randLimit_maxwell_multigrid < lb_loop_maxwell',isnan(bestParams),isinf(bestParams)],2));
                                                   lb_randLimit_maxwell_multigrid(badInds) = lb_loop_maxwell(badInds)';
                                                end

                                                if strcmp(elasticSetting_maxwell,'y')
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                        modulusInds = horzcat(1,(2:2:length(bestParams)-1));
                                                    else
                                                        modulusInds = horzcat(1,(2:2:length(bestParams)));
                                                    end
                                                else
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                        modulusInds = (1:2:length(bestParams)-1);
                                                    else
                                                        modulusInds = (1:2:length(bestParams));
                                                    end
                                                end

                                                ub_randLimit_maxwell_multigrid(modulusInds) = log10(ub_randLimit_maxwell_multigrid(modulusInds));
                                                lb_randLimit_maxwell_multigrid(modulusInds) = log10(lb_randLimit_maxwell_multigrid(modulusInds));

                                                resnorm_dist_temp = resnorm_dist_schapery_area;
                                                resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                                                [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                                                bestParams = beta_dist_schapery_area(:,idx);
                                                stdBand = std(beta_dist_schapery_area,0,2);
                                                ub_randLimit_schapery_multigrid = (bestParams + stdBand);
                                                lb_randLimit_schapery_multigrid = (bestParams - stdBand);

                                                if any(ub_randLimit_schapery_multigrid > ub_loop_maxwell','all') || any(lb_randLimit_schapery_multigrid < lb_loop_maxwell','all') || any(isnan(bestParams),'all') || any(isinf(bestParams),'all')
                                                   badInds = logical(sum([ub_randLimit_schapery_multigrid > ub_loop_maxwell',isnan(bestParams),isinf(bestParams)],2));
                                                   ub_randLimit_schapery_multigrid(badInds) = ub_loop_maxwell(badInds)';
                                                   badInds = logical(sum([lb_randLimit_schapery_multigrid < lb_loop_maxwell',isnan(bestParams),isinf(bestParams)],2));
                                                   lb_randLimit_schapery_multigrid(badInds) = lb_loop_maxwell(badInds)';
                                                end

                                                ub_randLimit_schapery_multigrid(modulusInds) = log10(ub_randLimit_schapery_multigrid(modulusInds));
                                                lb_randLimit_schapery_multigrid(modulusInds) = log10(lb_randLimit_schapery_multigrid(modulusInds));

                                            end
                                        else
                                            multiGridFlag = 0;
                                        end
                                    else
                                        multiGridFlag = 0;
                                    end
                                    multiGrid_iters = multiGrid_iters + 1;

                                end % End MultiGrid While

                            case 'nls'

                                % Note: a variable named "manualParams" must exist for
                                % this setting to be turned on (given a value of 1).
                                manualVoigtInputParam = 0;
                                manualMaxwellInputParam = 0;

                                if weightedFit
                                    % New weighting based on dataset
%                                     weightVectorVoigt = (1./h_norm_log_loop);
%                                     weightVectorMaxwell = (1./F_norm_log_loop);
%                                     weightVectorSchapery = (1./F_log_loop);
%                                     weightVectorDataVoigt = (1./h_norm_log_loop).*(1./(1-exp(t_plot_loop.*((nu_sample-1)./(dt)))));
%                                     weightVectorDataMaxwell = (1./F_norm_log_loop).*((1-exp(t_plot_loop.*((nu_sample-1)./(dt)))));
%                                     weightVectorDataSchapery = (1./F_log_loop).*(1+nu_sample).*((1-exp(t_plot_loop.*((nu_sample-1)./(dt)))));
            
                                    % time weighting
                                    timeweight = (t_plot_loop./(max(t_plot_loop)*0.9));
                                    timeweight(timeweight > 1) = 1;
                                    weightVectorVoigt = timeweight;
                                    weightVectorMaxwell = timeweight;
                                    weightVectorSchapery = timeweight;
                                    weightVectorDataVoigt = timeweight;
                                    weightVectorDataMaxwell = timeweight;
                                    weightVectorDataSchapery = timeweight;
            
                                else
                                    weightVectorVoigt = 1;
                                    weightVectorMaxwell = 1;
                                    weightVectorSchapery = 1;
                                    weightVectorDataVoigt = 1;
                                    weightVectorDataMaxwell = 1;
                                    weightVectorDataSchapery = 1;
                                end

            %                     %%%%%%%%%%%%%%%%%%%%%%%
            %                     % The result using "optimal" parameters.
            %                     % You need to convert Maxwell to Voigt parameters if
            %                     % your simulation data uses the Maxwell model. The same
            %                     % is true in reverse.
            %                     if elasticSetting_maxwell == 'y'
            %                         elParam = 1;
            %                     else
            %                         elParam = 0;
            %                     end
            % %                     vPlotParams = convertParams(manualParams,'voigt',elasticSetting_maxwell,fluidSetting_maxwell,elasticSetting,fluidSetting,(length(manualParams)-elParam)/2,x_fit_loop(:,1),dt,st,minTimescale);
            %                     vPlotParams = convertParamsCollocation(manualParams,'voigt',elasticSetting_maxwell,fluidSetting_maxwell,elasticSetting,fluidSetting,(length(manualParams)-elParam)/2,padSize,x_fit_loop(:,1),dt,st,minTimescale,advancedLog,forwardFitTimescale)';
            % 
            %                     figure
            %                     scatter(t_log_loop,weightVectorDataVoigt.*y_fit_loop,'ro')
            %                     hold on
            %                     scatter(t_log_loop,weightVectorDataMaxwell.*y_fit_maxwell_loop,'ko')
            %                     scatter(t_log_loop,weightVectorDataSchapery.*y_fit_area_loop,'bo')
            %                     
            %                     plot(t_log_loop,weightVectorVoigt.*F_conv_wrapper_loop(vPlotParams,x_fit_loop),'r')
            %                     plot(t_log_loop,0.1*weightedFit+weightVectorMaxwell.*h_conv_maxwell_wrapper_loop(manualParams,maxwellInputs_loop),'k')
            %                     plot(t_log_loop,weightVectorSchapery.*h_conv_schapery_area_wrapper_loop(manualParams,x_fit_area_loop),'b')
            %                     set(gca,'yscale','log')
            %                     set(gca,'xscale','log')
            %                     hold off
            %                     %%%%%%%%%%%%%%%%%%%%%%%

                                if fitSecondDatastream
                                    if strcmp(smoothData, 'n')
                                        inputsForce = [dataStruct(indShift+i_loop).t_r,...
                                            dataStruct(indShift+i_loop).h_r.^1.5];
                                        inputsForceSchapery = [dataStruct(indShift+i_loop).t_r,...
                                            sqrt(2*r_tip*dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>0)...
                                            - (dataStruct(indShift+i_loop).h_r(dataStruct(indShift+i_loop).t_r>0).^2))./r_tip];
                                    else
                                        inputsForce = [dataStruct(indShift+i_loop).t_r_smooth,...
                                            dataStruct(indShift+i_loop).h_r_smooth.^1.5];
                                        inputsForceSchapery = [dataStruct(indShift+i_loop).t_r_smooth(),...
                                            sqrt(2*r_tip*dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>0)...
                                            - (dataStruct(indShift+i_loop).h_r_smooth(dataStruct(indShift+i_loop).t_r_smooth>0).^2))./r_tip];
                                    end

                                    alpha = (8*sqrt(r_tip))/(3*(1-nu_sample));

                                    if weightedFit
                                        % Use the inverse of the data variance
                                        weightVector_alt = 1./var(dataStruct(indShift+i_loop).F_r,1);
                                    else
                                        weightVector_alt = 1;
                                    end

                                    % Make second data stream functions
                                    voigtForce = @(c) log_scale(alpha.*(h_conv_maxwell_loop(c,inputsForce)),timeVec,dt,st);
                                    maxwellForce = @(c) log_scale(alpha.*(h_conv_maxwell_loop(c,inputsForce)),timeVec,dt,st);
                                    schaperyMaxwellForce = @(c) log_scale(h_conv_schapery_area_loop(c,inputsForce),timeVec,dt,st);

                                    % Make wrapper to stack both the integral and
                                    % second data stream functions
                                    voigtFitWrapper = @(c,inputs) vertcat(weightVectorVoigt.*F_conv_wrapper_loop(c,inputs), weightVector_alt.*voigtForce(convertParams(c,'maxwell',elasticSetting,fluidSetting,elasticSetting_maxwell,fluidSetting_maxwell,i,maxwellInputs_loop(:,1),dt,st,minTimescale)));
                                    maxwellFitWrapper = @(c,inputs) vertcat(weightVectorMaxwell.*h_conv_maxwell_wrapper_loop(c,inputs), weightVector_alt.*maxwellForce(c));
                                    schaperyFitWrapper = @(c,inputs) vertcat(weightVectorSchapery.*h_conv_schapery_area_wrapper_loop(c,inputs), weightVector_alt.*schaperyMaxwellForce(c));

                                    % New ydata series
                                    y_fit_loop_alt = vertcat(weightVectorVoigt.*y_fit_loop, log_scale(weightVector_alt.*dataStruct(indShift+i_loop).F_r,timeVec,dt,st));
                                    y_fit_maxwell_loop_alt = vertcat(weightVectorMaxwell.*y_fit_maxwell_loop, log_scale(weightVector_alt.*dataStruct(indShift+i_loop).F_r,timeVec,dt,st));
                                    y_fit_area_loop_alt = vertcat(weightVectorSchapery.*y_fit_area_loop, log_scale(weightVector_alt.*dataStruct(indShift+i_loop).F_r,timeVec,dt,st));
                                end

                                ydata_loop = F_r_loop;
                                ydata_voigt_loop = h_r_loop;
                                timeVec_loop = x_fit_loop(:,1);
                                inputs_loop = [x_fit_loop(:,1),...
                                    F_r_loop];
                                inputsForce = [x_fit_loop(:,1),...
                                    h_r_loop.^1.5];
                                inputsForce_schapery = [x_fit_loop(x_fit_loop(:,1)>0,1),...
                                    sqrt(2*r_tip*h_r_loop(x_fit_loop(:,1)>0)...
                                    - (h_r_loop(x_fit_loop(:,1)>0).^2))./r_tip];
                                alphaModel = (8*sqrt(r_tip))/(3*(1-nu_sample));

                                % Parallel Loop for Random Search
                                parfor j = 1:loopLim
                                    warning('off');
                                    lastwarn(''); % Clear last warning

                                    try
                                        % Initialize
                                        newInds_linear = [];
                                        newInds_conv = [];
                                        newInds_maxwell = [];
                                        newInds_schapery = [];
                                        beta0_linear_loop = [];
                                        beta0_conv_loop = [];
                                        beta0_maxwell_loop = [];
                                        beta0_schapery_loop = [];
                                        ub_linear_loop = [];
                                        lb_linear_loop = [];
                                        ub_conv_loop = [];
                                        lb_conv_loop = [];
                                        ub_maxwell_loop = [];
                                        lb_maxwell_loop = [];
                                        ub_schapery_loop = [];
                                        lb_schapery_loop = [];
                                        ub_randLimit = [];
                                        lb_randLimit = [];
                                        ub_randLimit_maxwell = [];
                                        lb_randLimit_maxwell = [];
                                        ub_randLimit_schapery = [];
                                        lb_randLimit_schapery = [];
                                        beta0_temp_loop = [];
                                        limScale = [];
                                        tauInds = [];
                                        oldTau = [];
                                        modulusInds = [];
                                        oldModuli = [];
                                        weightArrayConv = [];
                                        weightArrayMaxwell = [];
                                        weightArraySchapery = [];

                                        % Find the new indices
                                        newInds_linear = isnan(beta0_linear);
                                        newInds_conv = isnan(beta0_conv);
                                        newInds_maxwell = isnan(beta0_maxwell);
                                        newInds_schapery = isnan(beta0_schapery);

                                        % Grab the base values for the beta0 cases
                                        beta0_linear_loop = beta0_linear;
                                        beta0_conv_loop = beta0_conv;
                                        beta0_maxwell_loop = beta0_maxwell;
                                        beta0_schapery_loop = beta0_schapery;

                                        % Control the bounds AND the initial parameters
                                        ub_randLimit = ub_randLimit_multigrid;
                                        lb_randLimit = lb_randLimit_multigrid;
                                        ub_randLimit_maxwell = ub_randLimit_maxwell_multigrid;
                                        lb_randLimit_maxwell = lb_randLimit_maxwell_multigrid;
                                        ub_randLimit_schapery = ub_randLimit_schapery_multigrid;
                                        lb_randLimit_schapery = lb_randLimit_schapery_multigrid;
                                        ub_fluidityRandLimit_init = 1e-20;
                                        limScale = getfield(logspace(5,0,n_terms),{i})*relaxationFactor;

                                        % Include or remove elastic term from bound
                                        % loosening that happens below. Use 1 for
                                        % inclusion, empty brackets [] for ignoring.
                                        elasticInd = 1;

            %                             oldModReductionFactor = limScale*0.5;
                                        oldModReductionFactor = 1;
            %                             oldModReductionFactor = 1/relaxationFactor;

                                        % Pick the right functions to use
                                        if fitLog
                                            if ~fitWithLSQ
                                                F_conv_func = @(c,inputs) F_conv_wrapper_loop(c,inputs);
                                                h_maxwell_func = @(c,inputs) h_conv_maxwell_wrapper_loop(c,inputs);
                                                h_schapery_func = @(c,inputs) h_conv_schapery_area_wrapper_loop(c,inputs);
                                            else
                                                if ~useJacobian
                                                    F_conv_func = @(c,inputs) F_conv_wrapper_loop(c,inputs);
                                                    h_maxwell_func = @(c,inputs) h_conv_maxwell_wrapper_loop(c,inputs);
                                                    h_schapery_func = @(c,inputs) h_conv_schapery_area_wrapper_loop(c,inputs);
                                                else
                                                    F_conv_func = @(c,inputs) F_conv_loop(c,inputs);
                                                    h_maxwell_func = @(c,inputs) h_conv_maxwell_loop(c,inputs);
                                                    h_schapery_func = @(c,inputs) h_conv_schapery_area_loop(c,inputs);
                                                end
                                            end
                                        else
                                            if ~fitWithLSQ
                                                F_conv_func = @(c,inputs) F_conv_loop(c,inputs);
                                                h_maxwell_func = @(c,inputs) h_conv_maxwell_loop(c,inputs);
                                                h_schapery_func = @(c,inputs) h_conv_schapery_area_loop(c,inputs);
                                            else
                                                if ~useJacobian
                                                    F_conv_func = @(c,inputs) F_conv_loop(c,inputs);
                                                    h_maxwell_func = @(c,inputs) h_conv_maxwell_loop(c,inputs);
                                                    h_schapery_func = @(c,inputs) h_conv_schapery_area_loop(c,inputs);
                                                else
                                                    F_conv_func = @(c,inputs) F_conv_loop(c,inputs);
                                                    h_maxwell_func = @(c,inputs) h_conv_maxwell_loop(c,inputs);
                                                    h_schapery_func = @(c,inputs) h_conv_schapery_area_loop(c,inputs);
                                                end
                                            end
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

                                        opts = optimoptions('fmincon','Algorithm','active-set',...
                                            'MaxFunctionEvaluations', n_maxIterations,...
                                            'MaxIterations', n_maxIterations,...
                                            'Display','none',...
                                            'FiniteDifferenceType','central',...
                                            'FunctionTolerance', 0,...
                                            'OptimalityTolerance', 0,...
                                            'StepTolerance', 0);

                                        if ~bruteForceBeta

                                            % Create random starting point
                                            [beta0_linear_loop,tauInds,modulusInds] = makeRandomParams(beta0_linear,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,newInds_linear);
                                            ub_linear_loop = ub_loop;
                                            lb_linear_loop = lb_loop;

                                            if ~manualVoigtInputParam

                                                % Create random starting point
                                                [beta0_conv_loop,tauInds,modulusInds] = makeRandomParams(beta0_conv,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,newInds_conv);

                                                weightArrayConv = ones(size(beta0_conv_loop));
                                                if scaleParams
                                                    numTau = length(tauInds);
                                                    if strcmp(fluidSetting,'y')
                                                       numTau = numTau-1;
                                                    end
                                                    weightArrayConv(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    if i == 1
                                                        weightArrayConv(modulusInds) = 10.^((max(ub_randLimit(modulusInds))-min(lb_randLimit(modulusInds)))/2);
                                                    else
                                                        switch lower(paramScale)
                                                            case 'equal'
                                                                weightArrayConv(modulusInds) = 10.^(floor(log10(mean(beta0_conv(~newInds_conv)))));
                                                            case 'individual'
                                                                weightArrayConv(intersect(find(~newInds_conv),modulusInds)) = 10.^(floor(log10(beta0_conv(intersect(find(~newInds_conv),modulusInds)))));
                                                        end
                                                    end
                                                end

                                            else
                                                % For testing against ELopezSim
                                                beta0_conv_loop = beta0_conv;
                                                if strcmp(elasticSetting,'y')
                                                    if strcmp(fluidSetting,'y')
                                                        beta0_conv_loop(1) = manualParams(1);
                                                        beta0_conv_loop(2:2:length(beta0_conv)-1) = manualParams(2:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(3:2:length(beta0_conv)-1) = manualParams(3:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(length(beta0_conv)) = manualParams(length(beta0_conv));
                                                    else
                                                        beta0_conv_loop(1) = manualParams(1);
                                                        beta0_conv_loop(2:2:length(beta0_conv)) = manualParams(2:2:length(beta0_conv));
                                                        beta0_conv_loop(3:2:length(beta0_conv)) = manualParams(3:2:length(beta0_conv));
                                                    end
                                                else
                                                    if strcmp(fluidSetting,'y')
                                                        beta0_conv_loop(1:2:length(beta0_conv)-1) = manualParams(1:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(2:2:length(beta0_conv)-1) = manualParams(2:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(length(beta0_conv)) = manualParams(length(beta0_conv));
                                                    else
                                                        beta0_conv_loop(1:2:length(beta0_conv)) = manualParams(1:2:length(beta0_conv));
                                                        beta0_conv_loop(2:2:length(beta0_conv)) = manualParams(2:2:length(beta0_conv));
                                                    end
                                                end

                                                weightArrayConv = ones(size(beta0_conv_loop));
                                            end

                                            % Write weighted function
                                            F_conv_loop_weight = weightParams(F_conv_func,weightArrayConv);
                                            ub_conv_loop = ub_loop;
                                            lb_conv_loop = lb_loop;
                                                
                                            if ~manualMaxwellInputParam

                                                % Create random starting point
                                                [beta0_maxwell_loop,tauInds,modulusInds] = makeRandomParams(beta0_maxwell,ub_randLimit_maxwell,lb_randLimit_maxwell,elasticSetting_maxwell,fluidSetting_maxwell,newInds_maxwell);

                                                weightArrayMaxwell = ones(size(beta0_maxwell_loop));
                                                if scaleParams
                                                    numTau = length(tauInds);
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                       numTau = numTau-1;
                                                    end
                                                    weightArrayMaxwell(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    if i == 1
                                                        weightArrayMaxwell(modulusInds) = 10.^((max(ub_randLimit_maxwell(modulusInds))-min(lb_randLimit_maxwell(modulusInds)))/2);
                                                    else
                                                        switch lower(paramScale)
                                                            case 'equal'
                                                                weightArrayMaxwell(modulusInds) = 10.^(floor(log10(mean(beta0_maxwell(~newInds_maxwell)))));
                                                            case 'individual'
                                                                weightArrayMaxwell(intersect(find(~newInds_conv),modulusInds)) = 10.^(floor(log10(beta0_conv(intersect(find(~newInds_conv),modulusInds)))));
                                                        end
                                                    end
                                                end

                                            else
                                                % For testing against ELopezSim
                                                beta0_maxwell_loop = beta0_maxwell;
                                                if strcmp(elasticSetting_maxwell,'y')
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                        beta0_maxwell_loop(1) = manualParams(1);
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)-1) = manualParams(2:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(3:2:length(beta0_maxwell)-1) = manualParams(3:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(length(beta0_maxwell)) = manualParams(length(beta0_maxwell));
                                                    else
                                                        beta0_maxwell_loop(1) = manualParams(1);
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)) = manualParams(2:2:length(beta0_maxwell));
                                                        beta0_maxwell_loop(3:2:length(beta0_maxwell)) = manualParams(3:2:length(beta0_maxwell));
                                                    end
                                                else
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                        beta0_maxwell_loop(1:2:length(beta0_maxwell)-1) = manualParams(1:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)-1) = manualParams(2:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(length(beta0_maxwell)) = manualParams(length(beta0_maxwell));
                                                    else
                                                        beta0_maxwell_loop(1:2:length(beta0_maxwell)) = manualParams(1:2:length(beta0_maxwell));
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)) = manualParams(2:2:length(beta0_maxwell));
                                                    end
                                                end

                                                weightArrayMaxwell = ones(size(beta0_maxwell_loop));
                                            end

                                            % Write weighted function
                                            h_conv_maxwell_loop_weight = weightParams(h_maxwell_func,weightArrayMaxwell);
                                            ub_maxwell_loop = ub_loop_maxwell;
                                            lb_maxwell_loop = lb_loop_maxwell;

                                            % Create random starting point
                                            [beta0_schapery_loop,tauInds,modulusInds] = makeRandomParams(beta0_schapery,ub_randLimit_schapery,lb_randLimit_schapery,elasticSetting_maxwell,fluidSetting_maxwell,newInds_schapery);

                                            weightArraySchapery = ones(size(beta0_schapery_loop));
                                            if scaleParams
                                                numTau = length(tauInds);
                                                if strcmp(fluidSetting_maxwell,'y')
                                                   numTau = numTau-1;
                                                end
                                                weightArraySchapery(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                if i == 1
                                                    weightArraySchapery(modulusInds) = 10.^((max(ub_randLimit_schapery(modulusInds))-min(lb_randLimit_schapery(modulusInds)))/2);
                                                else
                                                    switch lower(paramScale)
                                                        case 'equal'
                                                            weightArraySchapery(modulusInds) = 10.^(floor(log10(mean(beta0_maxwell(~newInds_schapery)))));
                                                        case 'individual'
                                                            weightArraySchapery(intersect(find(~newInds_conv),modulusInds)) = 10.^(floor(log10(beta0_conv(intersect(find(~newInds_conv),modulusInds)))));
                                                    end
                                                end
                                            end

                                            % Write weighted function
                                            h_conv_schapery_area_loop_weight = weightParams(h_schapery_func,weightArraySchapery);
                                            ub_schapery_loop = ub_loop_maxwell;
                                            lb_schapery_loop = lb_loop_maxwell;

                                        else

                                            if ~manualVoigtInputParam
                                                beta0_conv_loop = betaGrid(j,:);
                                            else
                                                beta0_conv_loop = beta0_conv;
                                                if strcmp(elasticSetting,'y')
                                                    if strcmp(fluidSetting,'y')
                                                        beta0_conv_loop(1) = manualParams(1);
                                                        beta0_conv_loop(2:2:length(beta0_conv)-1) = manualParams(2:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(3:2:length(beta0_conv)-1) = manualParams(3:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(length(beta0_conv)) = manualParams(length(beta0_conv));
                                                    else
                                                        beta0_conv_loop(1) = manualParams(1);
                                                        beta0_conv_loop(2:2:length(beta0_conv)) = manualParams(2:2:length(beta0_conv));
                                                        beta0_conv_loop(3:2:length(beta0_conv)) = manualParams(3:2:length(beta0_conv));
                                                    end
                                                else
                                                    if strcmp(fluidSetting,'y')
                                                        beta0_conv_loop(1:2:length(beta0_conv)-1) = manualParams(1:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(2:2:length(beta0_conv)-1) = manualParams(2:2:length(beta0_conv)-1);
                                                        beta0_conv_loop(length(beta0_conv)) = manualParams(length(beta0_conv));
                                                    else
                                                        beta0_conv_loop(1:2:length(beta0_conv)) = manualParams(1:2:length(beta0_conv));
                                                        beta0_conv_loop(2:2:length(beta0_conv)) = manualParams(2:2:length(beta0_conv));
                                                    end
                                                end
                                            end
                                            ub_conv_loop = ub_loop;
                                            lb_conv_loop = lb_loop;

                                            if ~manualMaxwellInputParam
                                                beta0_maxwell_loop = betaGrid_maxwell(j,:);
                                            else
                                                beta0_maxwell_loop = beta0_maxwell;
                                                if strcmp(elasticSetting_maxwell,'y')
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                        beta0_maxwell_loop(1) = manualParams(1);
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)-1) = manualParams(2:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(3:2:length(beta0_maxwell)-1) = manualParams(3:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(length(beta0_maxwell)) = manualParams(length(beta0_maxwell));
                                                    else
                                                        beta0_maxwell_loop(1) = manualParams(1);
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)) = manualParams(2:2:length(beta0_maxwell));
                                                        beta0_maxwell_loop(3:2:length(beta0_maxwell)) = manualParams(3:2:length(beta0_maxwell));
                                                    end
                                                else
                                                    if strcmp(fluidSetting_maxwell,'y')
                                                        beta0_maxwell_loop(1:2:length(beta0_maxwell)-1) = manualParams(1:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)-1) = manualParams(2:2:length(beta0_maxwell)-1);
                                                        beta0_maxwell_loop(length(beta0_maxwell)) = manualParams(length(beta0_maxwell));
                                                    else
                                                        beta0_maxwell_loop(1:2:length(beta0_maxwell)) = manualParams(1:2:length(beta0_maxwell));
                                                        beta0_maxwell_loop(2:2:length(beta0_maxwell)) = manualParams(2:2:length(beta0_maxwell));
                                                    end
                                                end
                                            end
                                            ub_maxwell_loop = ub_loop_maxwell;
                                            lb_maxwell_loop = lb_loop_maxwell;

                                            beta0_schapery_loop = betaGrid_schapery(j,:);
                                            ub_schapery_loop = ub_loop_maxwell;
                                            lb_schapery_loop = lb_loop_maxwell;

                                            beta0_linear_loop = betaGrid(j,:);
                                            ub_linear_loop = ub_loop;
                                            lb_linear_loop = lb_loop;

                                            % Write weighted functions
                                            weightArrayConv = ones(size(beta0_conv_loop));
                                            F_conv_loop_weight = weightParams(F_conv_func,weightArrayConv);
                                            weightArrayMaxwell = ones(size(beta0_maxwell_loop));
                                            h_conv_maxwell_loop_weight = weightParams(h_maxwell_func,weightArrayMaxwell);
                                            weightArraySchapery = ones(size(beta0_schapery_loop));
                                            h_conv_schapery_area_loop_weight = weightParams(h_schapery_func,weightArraySchapery);

                                        end

                                        % One last check of the parameters to try...
                                        if any(isnan(beta0_conv_loop))
                                            tempInd = isnan(beta0_conv_loop);
                                            beta0_conv_loop(tempInd) = lb_conv_loop(tempInd)+rand.*(ub_conv_loop(tempInd)-lb_conv_loop(tempInd));
                                        end
                                        if any(isnan(beta0_maxwell_loop))
                                            tempInd = isnan(beta0_maxwell_loop);
                                            beta0_maxwell_loop(tempInd) = lb_maxwell_loop(tempInd)+rand.*(ub_maxwell_loop(tempInd)-lb_maxwell_loop(tempInd));
                                        end
                                        if any(isnan(beta0_schapery_loop))
                                            tempInd = isnan(beta0_schapery_loop);
                                            beta0_schapery_loop(tempInd) = lb_schapery_loop(tempInd)+rand.*(ub_schapery_loop(tempInd)-lb_schapery_loop(tempInd));
                                        end

            %                             if any(isnan([beta0_conv_loop;beta0_maxwell_loop;beta0_schapery_loop]))
            %                                disp('NaN in the beta0 arrays!'); 
            %                             end

                                        if timeOperation
                                            Par.tic;
                                        end

                                        if ~fitWithLSQ

                                            % Fit the GKV model
            %                                 Fsumsquares = @(c) sum(((F_conv_wrapper_loop(c,x_fit_loop) - abs(y_fit))./abs(y_fit_loop)).^2);
            %                                 Fsumsquares = @(c) sum((F_conv_wrapper_loop(c,x_fit_loop) - y_fit_loop).^2);
                                            Fsumsquares = @(c) sum((weightVectorVoigt.*F_conv_loop_weight(c,x_fit_loop) - weightVectorDataVoigt.*y_fit_loop).^2);

                                            opts.StepTolerance = scaling;
            %                                 opts.Display = 'none';
            %                                 opts.OutputFcn = @outfun;

                                            [betatemp,ressquared,eflag,outputu] = ...
                                                fmincon(Fsumsquares,beta0_conv_loop,[],[],[],[],lb_conv_loop,ub_conv_loop,[],opts);
                                            betatemp = weightArrayConv.*betatemp;

                                            opts.StepTolerance = 0;
            %                                 opts.Display = 'none';
            %                                 opts.OutputFcn = [];

                                            residual = (F_conv_func(betatemp,x_fit_loop) - y_fit_loop);

                                        else

                                            if ~useJacobian
                                                lsqoptions.StepTolerance = scaling;
                %                                 lsqoptions.Display = 'none';
                %                                 lsqoptions.MaxFunctionEvaluations = 1e3;
                %                                 lsqoptions.OutputFcn = @outfun;

                                                if ~fitSecondDatastream
                                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(@(c,inputs) weightVectorVoigt.*F_conv_loop_weight(c,inputs),beta0_conv_loop,x_fit_loop,weightVectorDataVoigt.*y_fit_loop,lb_conv_loop,ub_conv_loop,lsqoptions);
                                                    betatemp = weightArrayConv.*betatemp;
                                                else
                                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(voigtFitWrapper,beta0_conv_loop,x_fit_loop,y_fit_loop_alt,lb_conv_loop,ub_conv_loop,lsqoptions);
                                                end
                                                lsqoptions.StepTolerance = 0;
                %                                 lsqoptions.Display = 'none';
                %                                 lsqoptions.MaxFunctionEvaluations = n_maxIterations;
                %                                 lsqoptions.OutputFcn = [];
                                            else
                                                lsqoptions.SpecifyObjectiveGradient = true;

                                                integralCalc = @(c,inputs,dQdI) subref_loop(convnfft(dQdI(c,inputs), inputs(:,2),'full'))*dt;
                                                [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(@(c,inputs) F_conv_wrapper_Jacobian(c,inputs,F_conv_loop_weight,integralCalc,elasticSetting,fluidSetting,dt,st,weightArrayConv,symbolicJacobian,F_norm_loop), beta0_conv_loop, x_fit_loop, y_fit_loop, lb_conv_loop, ub_conv_loop, lsqoptions);
                                                betatemp = weightArrayConv.*betatemp;

                                                lsqoptions.SpecifyObjectiveGradient = false;
                                            end

                                        end

                                        if timeOperation
                                            convTime(j) = Par.toc;
                                        end

                                        if ~ignoreWarns
                                            [warnMsg, warnId] = lastwarn;
                                            if ~isempty(warnMsg)
                                                if contains(warnMsg,'inaccurate')
                                                    beta_dist(:,j) = NaN;
                                                    beta0_dist(:,j) = NaN;
                                                    resnorm_dist(j) = NaN;

                                                    beta_dist_schapery_area(:,j) = NaN;
                                                    beta0_dist_schapery_area(:,j) = NaN;
                                                    resnorm_dist_schapery_area(j) = NaN;

                                                    beta_dist_linear(:,j) = NaN;
                                                    beta0_dist_linear(:,j) = NaN;
                                                    resnorm_dist_linear(j) = NaN;

                                                    fprintf('Entered a NaN Row for ALL cases.\n');
                                                    continue;
                                                end
                                            end
                                        end

                                        beta_dist(:,j) = betatemp';
                                        beta0_dist(:,j) = beta0_conv_loop;
                                        if ~fitSecondDatastream 
            %                                 resnorm_dist(j) = sum((residual./y_fit_loop).^2);
            %                                 resnorm_dist(j) = sum((residual).^2);
                                            if forceResidual
                                                yModelIndentation = (1./alphaModel.*F_conv_loop(betatemp,inputs_loop)).^(2/3);
                                                resnorm_dist(j) = sum(((weightVectorVoigt.*yModelIndentation-weightVectorDataVoigt.*ydata_voigt_loop).^2)./movvar(ydata_voigt_loop,3))./(length(ydata_voigt_loop)-length(betatemp)); % Standard Error of Regression (S)
                                                voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',log_scale(yModelIndentation,x_fit_loop(:,1),dt,st)'];
                                            else
                                                resnorm_dist(j) = sum(((weightVectorVoigt.*F_conv_func(betatemp,x_fit_loop)-weightVectorDataVoigt.*y_fit_loop).^2)./movvar(y_fit_loop,3))./(length(y_fit_loop)-length(betatemp)); % Standard Error of Regression (S)
                                                voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',F_conv_wrapper_loop(betatemp,x_fit_loop)'];
                                            end
            %                                 resnorm_dist(j) = sum(abs(((F_conv_wrapper_loop(betatemp,x_fit_loop)-y_fit_loop))./y_fit_loop));
                                        else
                                            resnorm_dist(j) = sum((residual./y_fit_loop_alt).^2,'all');
                                            voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',F_conv_wrapper_loop(betatemp,x_fit_loop)'];
                                        end

                                        if timeOperation
                                            Par.tic;
                                        end

                                        if ~fitWithLSQ

                                            % Fit the Maxwell model
            %                                 Fsumsquares_maxwell = @(c) sum(abs(h_conv_maxwell_wrapper_loop(c,maxwellInputs_loop) - y_fit_maxwell_loop)./y_fit_maxwell_loop);
            %                                 Fsumsquares_maxwell = @(c) sum((h_conv_maxwell_wrapper_loop(c,maxwellInputs_loop) - y_fit_maxwell_loop).^2);
                                            Fsumsquares_maxwell = @(c) sum((weightVectorMaxwell.*h_conv_maxwell_loop_weight(c,maxwellInputs_loop) - weightVectorDataMaxwell.*y_fit_maxwell_loop).^2);

                                            opts.StepTolerance = scaling;
            %                                 opts.Display = 'none';
            %                                 opts.OutputFcn = @outfun;

                                            [betatemp,ressquared,eflag,outputu] = ...
                                                fmincon(Fsumsquares_maxwell,beta0_maxwell_loop,[],[],[],[],lb_maxwell_loop,ub_maxwell_loop,[],opts);
                                            betatemp = weightArrayMaxwell.*betatemp;

                                            opts.StepTolerance = 0;
            %                                 opts.Display = 'none';
            %                                 opts.OutputFcn = [];

                                            residual = (h_maxwell_func(betatemp,maxwellInputs_loop) - y_fit_maxwell_loop);

                                        else

                                            if ~useJacobian
                                                lsqoptions.StepTolerance = scaling;
                %                                 lsqoptions.Display = 'none';
                %                                 lsqoptions.MaxFunctionEvaluations = 1e3;
                %                                 lsqoptions.OutputFcn = @outfun;

                                                if ~fitSecondDatastream 
                %                                     [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                %                                         lsqcurvefit(@(c,inputs) 0.1*weightedFit+weightVectorMaxwell.*h_conv_maxwell_wrapper_loop(c,inputs),beta0_maxwell_loop,maxwellInputs_loop,weightVectorDataMaxwell.*y_fit_maxwell_loop,lb_maxwell_loop,ub_maxwell_loop,lsqoptions);
                                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(@(c,inputs) weightVectorMaxwell.*h_maxwell_func(c,inputs),beta0_maxwell_loop,maxwellInputs_loop,weightVectorDataMaxwell.*y_fit_maxwell_loop,lb_maxwell_loop,ub_maxwell_loop,lsqoptions);
                                                    betatemp = weightArrayMaxwell.*betatemp;
                                                else
                                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(maxwellFitWrapper,beta0_maxwell_loop,maxwellInputs_loop,y_fit_maxwell_loop_alt,lb_maxwell_loop,ub_maxwell_loop,lsqoptions);
                                                end

                                                if exitflag == 1
                                                    disp('Maxwell Model Converged!');   
                                                end

                                                lsqoptions.StepTolerance = 0;
                %                                 lsqoptions.Display = 'none';
                %                                 lsqoptions.MaxFunctionEvaluations = n_maxIterations;
                %                                 lsqoptions.OutputFcn = [];
                                            else
                                                lsqoptions.SpecifyObjectiveGradient = true;

                                                integralCalc = @(c,inputs,dQdI) subref_loop(convnfft(dQdI(c,inputs), inputs(:,2),'full'))*dt;
                                                [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(@(c,inputs) h_conv_wrapper_Jacobian(c,inputs,h_conv_maxwell_loop_weight,integralCalc,elasticSetting_maxwell,fluidSetting_maxwell,dt,st,weightArrayMaxwell,symbolicJacobian,h_norm_loop), beta0_maxwell_loop, maxwellInputs_loop, y_fit_maxwell_loop, lb_maxwell_loop, ub_maxwell_loop, lsqoptions);
                                                betatemp = weightArrayMaxwell.*betatemp;    

                                                lsqoptions.SpecifyObjectiveGradient = false;
                                            end

                                        end

                                        if timeOperation
                                            maxwellTime(j) = Par.toc;
                                        end

                                        if ~ignoreWarns
                                            [warnMsg, warnId] = lastwarn;
                                            if ~isempty(warnMsg)
                                                if contains(warnMsg,'inaccurate')
                                                    beta_dist(:,j) = NaN;
                                                    beta0_dist(:,j) = NaN;
                                                    resnorm_dist(j) = NaN;

                                                    fprintf('Entered a NaN Row.\n');
                                                    continue;
                                                end
                                            end
                                        end

                                        beta_dist_maxwell(:,j) = betatemp';
                                        beta0_dist_maxwell(:,j) = beta0_maxwell_loop;
                                        if ~fitSecondDatastream 
            %                                 resnorm_dist_maxwell(j) = sum((residual./y_fit_maxwell_loop).^2);
            %                                 resnorm_dist_maxwell(j) = sum((residual).^2);
                                            if forceResidual
                                                if strcmp(elasticSetting_maxwell,'y')
                                                    yModelForce = alphaModel.*(betatemp(1).*(inputsForce(:,2))+subref_loop(convnfft(G_func_maxwell_loop(betatemp,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                                else
                                                    yModelForce = alphaModel.*(subref_loop(convnfft(G_func_maxwell_loop(betatemp,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                                end
                                                resnorm_dist_maxwell(j) = sum(((weightVectorMaxwell.*yModelForce-weightVectorDataMaxwell.*ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(betatemp)); % Standard Error of Regression (S)
                                                maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',log_scale(yModelForce,maxwellInputs_loop(:,1),dt,st)'];
                                            else
                                                resnorm_dist_maxwell(j) = sum(((weightVectorMaxwell.*h_maxwell_func(betatemp,maxwellInputs_loop)-weightVectorDataMaxwell.*y_fit_maxwell_loop).^2)./movvar(y_fit_maxwell_loop,3))./(length(y_fit_maxwell_loop)-length(betatemp)); % Standard Error of Regression (S)
                                                maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',h_conv_maxwell_wrapper_loop(betatemp,maxwellInputs_loop)'];
                                            end
            %                                 resnorm_dist_maxwell(j) = sum(abs(((h_conv_maxwell_wrapper_loop(betatemp,maxwellInputs_loop)-y_fit_maxwell_loop))./y_fit_maxwell_loop));
                                        else
                                            resnorm_dist_maxwell(j) = sum((residual./y_fit_maxwell_loop_alt).^2,'all');
                                            maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',h_conv_maxwell_wrapper_loop(betatemp,maxwellInputs_loop)'];
                                        end

                                        if timeOperation
                                            Par.tic;
                                        end

                                        if ~fitWithLSQ

                                            % Fit the Maxwell model
            %                                 Fsumsquaresschapery = @(c) sum(abs(h_conv_schapery_area_wrapper_loop(c,x_fit_area_loop) - y_fit_area_loop)./y_fit_area_loop);
            %                                 Fsumsquaresschapery = @(c) sum((h_conv_schapery_area_wrapper_loop(c,x_fit_area_loop) - y_fit_area_loop).^2);
                                            Fsumsquaresschapery = @(c) sum((weightVectorSchapery.*h_conv_schapery_area_loop_weight(c,x_fit_area_loop) - weightVectorDataSchapery.*y_fit_area_loop).^2);

                                            opts.StepTolerance = scaling;
            %                                 opts.Display = 'none';
            %                                 opts.OutputFcn = @outfun;

                                            [betatemp,ressquared,eflag,outputu] = ...
                                                fmincon(Fsumsquaresschapery,beta0_schapery_loop,[],[],[],[],lb_schapery_loop,ub_schapery_loop,[],opts);
                                            betatemp = weightArraySchapery.*betatemp;

                                            opts.StepTolerance = 0;
            %                                 opts.Display = 'none';
            %                                 opts.OutputFcn = [];

                                            residual = (h_schapery_func(betatemp,x_fit_area_loop) - y_fit_area_loop);

                                        else

                                            if ~useJacobian
                                                lsqoptions.StepTolerance = scaling;
                %                                 lsqoptions.MaxFunctionEvaluations = 1e3;
                %                                 lsqoptions.Display = 'none';
                %                                 lsqoptions.OutputFcn = @outfun;

                                                if ~fitSecondDatastream
                                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(@(c,inputs) weightVectorSchapery.*h_schapery_func(c,inputs),beta0_schapery_loop,x_fit_area_loop,weightVectorDataSchapery.*y_fit_area_loop,lb_schapery_loop,ub_schapery_loop,lsqoptions);
                                                    betatemp = weightArraySchapery.*betatemp;
                                                else
                                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(schaperyFitWrapper,beta0_schapery_loop,x_fit_area_loop,y_fit_area_loop_alt,lb_schapery_loop,ub_schapery_loop,lsqoptions);
                                                end

                                                if exitflag == 1
                                                    disp('Schapery-Maxwell Model Converged!');
                                                end

                                                lsqoptions.StepTolerance = 0;
                %                                 lsqoptions.MaxFunctionEvaluations = n_maxIterations;
                %                                 lsqoptions.Display = 'none';
                %                                 lsqoptions.OutputFcn = @outfun;
                                            else
                                                lsqoptions.SpecifyObjectiveGradient = true;

                                                integralCalc = @(c,inputs,dQdI) subref_loop(convnfft(dQdI(c,inputs), inputs(:,2),'full'))*dt;
                                                [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                                        lsqcurvefit(@(c,inputs) h_conv_wrapper_Jacobian(c,inputs,h_conv_schapery_area_loop_weight,integralCalc,elasticSetting_maxwell,fluidSetting_maxwell,dt,st,weightArraySchapery,symbolicJacobian,h_norm_loop), beta0_schapery_loop, x_fit_area_loop, y_fit_area_loop, lb_schapery_loop, ub_schapery_loop, lsqoptions);
                                                betatemp = weightArraySchapery.*betatemp;

                                                lsqoptions.SpecifyObjectiveGradient = false;
                                            end

                                        end

                                        if timeOperation
                                            schaperyTime(j) = Par.toc;
                                        end

                                        if ~ignoreWarns
                                            [warnMsg, warnId] = lastwarn;
                                            if ~isempty(warnMsg)
                                                if contains(warnMsg,'inaccurate')
                                                    beta_dist_schapery_area(:,j) = NaN;
                                                    beta0_dist_schapery_area(:,j) = NaN;
                                                    resnorm_dist_schapery_area(j) = NaN;

                                                    beta_dist_linear(:,j) = NaN;
                                                    beta0_dist_linear(:,j) = NaN;
                                                    resnorm_dist_linear(j) = NaN;

                                                    fprintf('Entered a NaN Row for the Schapery and Linear cases.\n');
                                                    continue;
                                                end
                                            end
                                        end

                                        beta_dist_schapery_area(:,j) = betatemp';
                                        beta0_dist_schapery_area(:,j) = beta0_schapery_loop;
                                        if ~fitSecondDatastream
            %                                 resnorm_dist_schapery_area(j) = sum((residual./y_fit_area_loop).^2);
            %                                 resnorm_dist_schapery_area(j) = sum((residual).^2);
                                            if forceResidual
                                                yModelForce = h_conv_schapery_area_loop(betatemp,inputsForce_schapery);
                                                resnorm_dist_schapery_area(j) = sum(((weightVectorSchapery.*yModelForce-weightVectorDataSchapery.*ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(betatemp)); % Standard Error of Regression (S)
                                                schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',log_scale(yModelForce,x_fit_area_loop(:,1),dt,st)'];
                                            else
                                                resnorm_dist_schapery_area(j) = sum(((weightVectorSchapery.*h_schapery_func(betatemp,x_fit_area_loop)-weightVectorDataSchapery.*y_fit_area_loop).^2)./movvar(y_fit_area_loop,3))./(length(y_fit_area_loop)-length(betatemp)); % Standard Error of Regression (S)
                                                schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',h_conv_schapery_area_wrapper_loop(betatemp,x_fit_area_loop)'];
                                            end
            %                                 resnorm_dist_schapery_area(j) = sum(abs(((h_conv_schapery_area_wrapper_loop(betatemp,x_fit_area_loop)-y_fit_area_loop))./y_fit_area_loop));
                                        else
                                            resnorm_dist_schapery_area(j) = sum((residual./y_fit_area_loop_alt).^2,'all');
                                            schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',h_conv_schapery_area_wrapper_loop(betatemp,x_fit_area_loop)'];
                                        end

                                        % Fit the linear model
            %                             Fsumsquares_linear = @(c) sum(((F_linear_fit_loop(c,x_fit_loop) - abs(y_fit_linear))./abs(y_fit_linear)));
                                        Fsumsquares_linear = @(c) sum(abs(F_linear_fit_loop(c,x_fit_loop(:,1)) - y_fit_linear_loop));
                                        opts = optimoptions('fmincon','Algorithm','interior-point',...
                                            'Display','none','FiniteDifferenceType','central',...
                                            'FunctionTolerance', scaling,...
                                            'OptimalityTolerance',scaling,...
                                            'StepTolerance', scaling);
                                        [betatemp,ressquared,eflag,outputu] = ...
                                            fmincon(Fsumsquares_linear,beta0_linear_loop,[],[],[],[],lb_linear_loop,ub_linear_loop,[],opts);

                                        if ~ignoreWarns
                                            [warnMsg, warnId] = lastwarn;
                                            if ~isempty(warnMsg)
                                                if contains(warnMsg,'inaccurate')
                                                    beta_dist_linear(:,j) = NaN;
                                                    beta0_dist_linear(:,j) = NaN;
                                                    resnorm_dist_linear(j) = NaN;

                                                    fprintf('Entered a NaN Row for the Linear case.\n');
                                                    continue;
                                                end
                                            end
                                        end

                                        beta_dist_linear(:,j) = betatemp';
                                        beta0_dist_linear(:,j) = beta0_linear_loop;
                                        resnorm_dist_linear(j) = ressquared;

                                        if timeOperation
                                            % Track memory
                                            %ramused(j) = userMem.MemUsedMATLAB;
                                            ramused(j) = 0;
                                            %ramav(j) = systemMem.PhysicalMemory.Available;
                                            ramav(j) = 0;
                                            ramtime(j) = now;
                                        end

                                    catch ERROR

                                        beta_dist(:,j) = NaN;
                                        beta0_dist(:,j) = NaN;
                                        resnorm_dist(j) = NaN;

                                        beta_dist_maxwell(:,j) = NaN;
                                        beta0_dist_maxwell(:,j) = NaN;
                                        resnorm_dist_maxwell(j) = NaN;

                                        beta_dist_schapery_area(:,j) = NaN;
                                        beta0_dist_schapery_area(:,j) = NaN;
                                        resnorm_dist_schapery_area(j) = NaN;

                                        beta_dist_linear(:,j) = NaN;
                                        beta0_dist_linear(:,j) = NaN;
                                        resnorm_dist_linear(j) = NaN;

                                        if timeOperation
                                            convTime(j) = 1;
                                            maxwellTime(j) = 1;
                                            schaperyTime(j) = 1;
                                        end

                                        voigtData = [0,0];
                                        maxwellData = [0,0];
                                        schaperyData = [0,0];

%                                         fprintf('\nTerm Number %d Search: Entered a NaN Row for loop #%d because there was an ERROR:\n',i,j);
%                                         disp(ERROR.message)
                                    end

                                end

                            case {'gen-alg','pattern','particle-swarm','surrogate'}

                                printData = 0;
                                plotCost = 0;
                                plotFit = 0;
                                n_plot = 1000;

                                errorCalcVoigt = @(c,inputs) sum(((F_conv_wrapper_loop(c,inputs)-y_fit_loop).^2)./movvar(y_fit_loop,3))./(length(y_fit_loop)-length(c)); % Standard Error of Regression (S)
                                errorCalcMaxwell = @(c,inputs) sum(((h_conv_maxwell_wrapper_loop(c,inputs)-y_fit_maxwell_loop).^2)./movvar(y_fit_maxwell_loop,3))./(length(y_fit_maxwell_loop)-length(c)); % Standard Error of Regression (S)
                                errorCalcSchapery = @(c,inputs) sum(((h_conv_schapery_area_wrapper_loop(c,inputs)-y_fit_area_loop).^2)./movvar(y_fit_area_loop,3))./(length(y_fit_area_loop)-length(c)); % Standard Error of Regression (S)

                                if plotCost
                                    if exist('costPlot','var')
                                        clearvars costPlot
                                    end
                                    if exist('costPlot','var')
                                        figure(costPlot);
                                        clf;
                                    else
                                        costPlot = figure('position',[25 25 600 400]);
                                    end
                                    if exist('fitPlot','var')
                                        clearvars fitPlot
                                    end
                                    if exist('fitPlot','var')
                                        figure(fitPlot);
                                        clf;
                                    else
                                        fitPlot = figure('position',[125 25 600 400]);
                                    end
                                end

                                ydata_loop = F_r_loop;
                                ydata_voigt_loop = h_r_loop;
                                timeVec_loop = x_fit_loop(:,1);
                                inputs_loop = [x_fit_loop(:,1),...
                                    F_r_loop];
                                inputsForce = [x_fit_loop(:,1),...
                                    h_r_loop.^1.5];
                                inputsForce_schapery = [x_fit_loop(x_fit_loop(:,1)>0,1),...
                                    sqrt(2*r_tip*h_r_loop(x_fit_loop(:,1)>0)...
                                    - (h_r_loop(x_fit_loop(:,1)>0).^2))./r_tip];
                                alphaModel = (8*sqrt(r_tip))/(3*(1-nu_sample));

                                y_fit_GD_loop = h_norm_loop;
                                y_fit_maxwell_GD_loop = F_norm_loop;
                                y_fit_schapery_GD_loop = F_schapery_area_loop;

                                switch lower(optimizationMethod)
                                    case 'gen-alg'
                                        titleString = 'Genetic Algorithm';
                                    case 'pattern'
                                        titleString = 'Pattern Search';
                                    case 'particle-swarm'
                                        titleString = 'Particle Swarm';
                                    case 'surrogate'
                                        titleString = 'Surrogate';
                                end

                                parfor j = 1:loopLim
                                    warning('off');
                                    % Initialize
                                    limScale = [];
                                    tauInds = [];
                                    modulusInds = [];
                                    newParams_voigt = [];
                                    newParams_maxwell = [];
                                    newParams_schapery = [];
                                    newTau = [];
                                    newMod = [];

                                    % Find the new indices
                                    newInds_linear = isnan(beta0_linear);
                                    newInds_conv = isnan(beta0_conv);
                                    newInds_maxwell = isnan(beta0_maxwell);
                                    newInds_schapery = isnan(beta0_schapery);

                                    % Grab the base values for the beta0 cases
                                    beta0_linear_loop = beta0_linear;
                                    beta0_conv_loop = beta0_conv;
                                    beta0_maxwell_loop = beta0_maxwell;
                                    beta0_schapery_loop = beta0_schapery;

                                    % Control the bounds AND the initial parameters
                                    ub_randLimit = ub_randLimit_multigrid;
                                    lb_randLimit = lb_randLimit_multigrid;
                                    ub_randLimit_maxwell = ub_randLimit_maxwell_multigrid;
                                    lb_randLimit_maxwell = lb_randLimit_maxwell_multigrid;
                                    ub_randLimit_schapery = ub_randLimit_schapery_multigrid;
                                    lb_randLimit_schapery = lb_randLimit_schapery_multigrid;
                                    ub_fluidityRandLimit_init = 1e-20;
                                    limScale = getfield(logspace(5,0,n_terms),{i})*relaxationFactor;

                                    % Include or remove elastic term from bound
                                    % loosening that happens below. Use 1 for
                                    % inclusion, empty brackets [] for ignoring.
                                    elasticInd = 1;

        %                             oldModReductionFactor = limScale*0.5;
                                    oldModReductionFactor = 1;
        %                             oldModReductionFactor = 1/relaxationFactor;
        
                                    [tauInds,modulusInds] = getParamIndices(beta0_conv_loop,elasticSetting,fluidSetting);
                                    [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(beta0_maxwell_loop,elasticSetting_maxwell,fluidSetting_maxwell);
                                    
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

                                    if any(newInds_conv) && ~enforceGridGuess

                                        beta0_temp_loop = beta0_conv;
                                        if strcmp(elasticSetting,'y')
                                            if strcmp(fluidSetting,'y')
                                                beta0_temp_loop(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                                                beta0_temp_loop(2:2:end-1) = 10.^(((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1)))+lb_randLimit(2:2:end-1)));
                                                beta0_temp_loop(3:2:end-1) = ((ub_randLimit(3:2:end-1)-lb_randLimit(3:2:end-1)).*rand(size(ub_randLimit(3:2:end-1))) + lb_randLimit(3:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                                                beta0_temp_loop(2:2:end) = 10.^(((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end)))+lb_randLimit(2:2:end)));
                                                beta0_temp_loop(3:2:end) = ((ub_randLimit(3:2:end)-lb_randLimit(3:2:end)).*rand(size(ub_randLimit(3:2:end))) + lb_randLimit(3:2:end));
                                            end
                                        else
                                            if strcmp(fluidSetting,'y')
                                                beta0_temp_loop(1:2:end-1) = 10.^(((ub_randLimit(1:2:end-1)-lb_randLimit(1:2:end-1)).*rand(size(ub_randLimit(1:2:end-1)))+lb_randLimit(1:2:end-1)));
                                                beta0_temp_loop(2:2:end-1) = ((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1))) + lb_randLimit(2:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1:2:end) = 10.^(((ub_randLimit(1:2:end)-lb_randLimit(1:2:end)).*rand(size(ub_randLimit(1:2:end)))+lb_randLimit(1:2:end)));
                                                beta0_temp_loop(2:2:end) = ((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end))) + lb_randLimit(2:2:end));
                                            end
                                        end

                                        beta0_conv_loop(newInds_conv) = beta0_temp_loop(newInds_conv);
                                        ub_conv_loop = ub_loop;
                                        lb_conv_loop = lb_loop;

                                        % Enforce the original time
                                        % constant bounds
                                        [oldTau,~] = ismember(find(~newInds_conv),tauInds);
                                        if ~isempty(oldTau)
                                            ub_conv_loop(oldTau) = ub_loop(oldTau);
                                            lb_conv_loop(oldTau) = lb_loop(oldTau);
                                        end

                                        % Allow compliances to drop without restriction
                                        [oldMods,~] = ismember(find(~newInds_conv),modulusInds);

                                        beta0_conv_loop(oldMods) = beta0_conv_loop(oldMods).*oldModReductionFactor;

                                        % Enforce the original upper and lower
                                        % stiffness bounds
                                        if any(lb_conv_loop(oldMods) < lb_loop(oldMods))
                                            lb_conv_loop(lb_conv_loop < lb_loop) = lb_loop(lb_conv_loop < lb_loop);
                                        end
                                        if any(ub_conv_loop(oldMods) > ub_loop(oldMods))
                                            ub_conv_loop(ub_conv_loop > ub_loop) = ub_loop(ub_conv_loop > ub_loop);
                                        end

        %                                 ub_conv_loop(oldMods) = 1e100;
                                        ub_conv_loop(oldMods) = ub_loop(1);

                                        if any(beta0_conv_loop > ub_conv_loop) || any(beta0_conv_loop < lb_conv_loop)
                                           badInds = logical((beta0_conv_loop > ub_conv_loop) + (beta0_conv_loop < lb_conv_loop));
                                           beta0_conv_loop(badInds) = (ub_conv_loop(badInds)-lb_conv_loop(badInds)) .* rand( size(beta0_conv_loop(badInds)) ) + lb_conv_loop(badInds);
                                        end

                                    elseif any(newInds_conv) && enforceGridGuess

                                        ub_conv_loop = ub_loop;
                                        lb_conv_loop = lb_loop;

                                        % Enforce the original time
                                        % constant bounds
                                        [oldTau,~] = ismember(find(~newInds_conv),tauInds);
                                        if ~isempty(oldTau)
                                            ub_conv_loop(oldTau) = ub_loop(oldTau);
                                            lb_conv_loop(oldTau) = lb_loop(oldTau);
                                        end

                                        % Allow compliances to drop without restriction
                                        [oldMods,~] = ismember(find(~newInds_conv),modulusInds);

                                        newTau = tauInds(ismember(tauInds,find(newInds_conv)));
                                        newMod = modulusInds(ismember(modulusInds,find(newInds_conv)));

                                        tauStep = floor(sqrt(loopLim));
                                        modStep = floor((loopLim)^(1./(length(ub_loop(newMod)))))/(length(newTau)*tauStep);
                                        [betaGridJ,betaGridTau] = meshgrid(logspace(log10(lb_loop(newMod)),log10(ub_loop(newMod)),modStep),...
                                            logspace(log10(lb_loop(newTau)),log10(ub_loop(newTau)),tauStep));

                                        if j <= tauStep*modStep
                                            beta0_conv_loop(newInds_conv) = [betaGridJ(j);betaGridTau(j)];
                                        else
                                            beta0_conv_loop(newInds_conv) = [betaGridJ(end);betaGridTau(end)];
                                        end

                                        beta0_conv_loop(oldMods) = beta0_conv_loop(oldMods).*oldModReductionFactor;

                                        % Enforce the original upper and lower
                                        % stiffness bounds
                                        if any(lb_conv_loop(oldMods) < lb_loop(oldMods))
                                            lb_conv_loop(lb_conv_loop < lb_loop) = lb_loop(lb_conv_loop < lb_loop);
                                        end
                                        if any(ub_conv_loop(oldMods) > ub_loop(oldMods))
                                            ub_conv_loop(ub_conv_loop > ub_loop) = ub_loop(ub_conv_loop > ub_loop);
                                        end

        %                                 ub_conv_loop(oldMods) = 1e100;
                                        ub_conv_loop(oldMods) = ub_loop(1);

                                        if any(beta0_conv_loop > ub_conv_loop) || any(beta0_conv_loop < lb_conv_loop)
                                           badInds = logical((beta0_conv_loop > ub_conv_loop) + (beta0_conv_loop < lb_conv_loop));
                                           beta0_conv_loop(badInds) = (ub_conv_loop(badInds)-lb_conv_loop(badInds)) .* rand( size(beta0_conv_loop(badInds)) ) + lb_conv_loop(badInds);
                                        end

                                    else

                                        beta0_temp_loop = beta0_conv;
                                        if strcmp(elasticSetting,'y')
                                            if strcmp(fluidSetting,'y')
                                                beta0_temp_loop(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                                                beta0_temp_loop(2:2:end-1) = 10.^(((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1)))+lb_randLimit(2:2:end-1)));
                                                beta0_temp_loop(3:2:end-1) = ((ub_randLimit(3:2:end-1)-lb_randLimit(3:2:end-1)).*rand(size(ub_randLimit(3:2:end-1))) + lb_randLimit(3:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                                                beta0_temp_loop(2:2:end) = 10.^(((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end)))+lb_randLimit(2:2:end)));
                                                beta0_temp_loop(3:2:end) = ((ub_randLimit(3:2:end)-lb_randLimit(3:2:end)).*rand(size(ub_randLimit(3:2:end))) + lb_randLimit(3:2:end));
                                            end
                                        else
                                            if strcmp(fluidSetting,'y')
                                                beta0_temp_loop(1:2:end-1) = 10.^(((ub_randLimit(1:2:end-1)-lb_randLimit(1:2:end-1)).*rand(size(ub_randLimit(1:2:end-1)))+lb_randLimit(1:2:end-1)));
                                                beta0_temp_loop(2:2:end-1) = ((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1))) + lb_randLimit(2:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1:2:end) = 10.^(((ub_randLimit(1:2:end)-lb_randLimit(1:2:end)).*rand(size(ub_randLimit(1:2:end)))+lb_randLimit(1:2:end)));
                                                beta0_temp_loop(2:2:end) = ((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end))) + lb_randLimit(2:2:end));
                                            end
                                        end

                                        if strcmp(elasticSetting,'y') && fitElasticFirst
                                            beta0_temp_loop(1) = (bestElasticTerm)*oldModReductionFactor;
                                            ub_conv_loop = ub_loop;
                                            lb_conv_loop = lb_loop;
                                        else
                                            ub_conv_loop = ub_loop;
                                            lb_conv_loop = lb_loop;
                                        end

                                        beta0_conv_loop = beta0_temp_loop;

                                    end

                                    % Make new starting point (Maxwell)
                                    if any(newInds_maxwell) && ~enforceGridGuess
                                        beta0_temp_loop = beta0_maxwell;
                                        if strcmp(elasticSetting_maxwell,'y')
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_maxwell(1)-lb_randLimit_maxwell(1))*rand()+lb_randLimit_maxwell(1)));
                                                beta0_temp_loop(2:2:end-1) = 10.^(((ub_randLimit_maxwell(2:2:end-1)-lb_randLimit_maxwell(2:2:end-1)).*rand(size(ub_randLimit_maxwell(2:2:end-1)))+lb_randLimit_maxwell(2:2:end-1)));
                                                beta0_temp_loop(3:2:end-1) = ((ub_randLimit_maxwell(3:2:end-1)-lb_randLimit_maxwell(3:2:end-1)).*rand(size(ub_randLimit_maxwell(3:2:end-1))) + lb_randLimit_maxwell(3:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_maxwell(1)-lb_randLimit_maxwell(1))*rand()+lb_randLimit_maxwell(1)));
                                                beta0_temp_loop(2:2:end) = 10.^(((ub_randLimit_maxwell(2:2:end)-lb_randLimit_maxwell(2:2:end)).*rand(size(ub_randLimit_maxwell(2:2:end)))+lb_randLimit_maxwell(2:2:end)));
                                                beta0_temp_loop(3:2:end) = ((ub_randLimit_maxwell(3:2:end)-lb_randLimit_maxwell(3:2:end)).*rand(size(ub_randLimit_maxwell(3:2:end))) + lb_randLimit_maxwell(3:2:end));
                                            end
                                        else
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1:2:end-1) = 10.^(((ub_randLimit_maxwell(1:2:end-1)-lb_randLimit_maxwell(1:2:end-1)).*rand(size(ub_randLimit_maxwell(1:2:end-1)))+lb_randLimit_maxwell(1:2:end-1)));
                                                beta0_temp_loop(2:2:end-1) = ((ub_randLimit_maxwell(2:2:end-1)-lb_randLimit_maxwell(2:2:end-1)).*rand(size(ub_randLimit_maxwell(2:2:end-1))) + lb_randLimit_maxwell(2:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1:2:end) = 10.^(((ub_randLimit_maxwell(1:2:end)-lb_randLimit_maxwell(1:2:end)).*rand(size(ub_randLimit_maxwell(1:2:end)))+lb_randLimit_maxwell(1:2:end)));
                                                beta0_temp_loop(2:2:end) = ((ub_randLimit_maxwell(2:2:end)-lb_randLimit_maxwell(2:2:end)).*rand(size(ub_randLimit_maxwell(2:2:end))) + lb_randLimit_maxwell(2:2:end));
                                            end
                                        end

                                        beta0_maxwell_loop(newInds_maxwell) = beta0_temp_loop(newInds_maxwell);
                                        ub_maxwell_loop = ub_loop_maxwell;
                                        lb_maxwell_loop = lb_loop_maxwell;

                                        [oldMods,~] = ismember(find(~newInds_maxwell),modulusInds_maxwell);
                                        beta0_maxwell_loop(oldMods) = beta0_maxwell_loop(oldMods)./oldModReductionFactor;

                                        % Enforce the original upper and lower
                                        % stiffness bounds
                                        if any(lb_maxwell_loop(oldMods) < lb_loop_maxwell(oldMods))
                                            lb_maxwell_loop(lb_maxwell_loop < lb_loop_maxwell) = lb_loop_maxwell(lb_maxwell_loop < lb_loop_maxwell);
                                        end
                                        if any(ub_maxwell_loop(oldMods) > ub_loop_maxwell(oldMods))
                                            ub_maxwell_loop(ub_maxwell_loop > ub_loop_maxwell) = ub_loop_maxwell(ub_maxwell_loop > ub_loop_maxwell);
                                        end

                                        % Enforce the original time
                                        % constant bounds
                                        [oldTau,~] = ismember(find(~newInds_maxwell),tauInds_maxwell);
                                        if ~isempty(oldTau)
                                            ub_maxwell_loop(oldTau) = ub_loop_maxwell(oldTau);
                                            lb_maxwell_loop(oldTau) = lb_loop_maxwell(oldTau);
                                        end

                                        if any(beta0_maxwell_loop > ub_maxwell_loop) || any(beta0_maxwell_loop < lb_maxwell_loop)
                                           badInds = logical((beta0_maxwell_loop > ub_maxwell_loop) + (beta0_maxwell_loop < lb_maxwell_loop));
                                           beta0_maxwell_loop(badInds) = (ub_maxwell_loop(badInds)-lb_maxwell_loop(badInds)) .* rand( size(beta0_maxwell_loop(badInds)) ) + lb_maxwell_loop(badInds);
                                        end

                                    elseif any(newInds_maxwell) && enforceGridGuess

                                        ub_maxwell_loop = ub_loop_maxwell;
                                        lb_maxwell_loop = lb_loop_maxwell;

                                        % Enforce the original time
                                        % constant bounds
                                        [oldTau,~] = ismember(find(~newInds_maxwell),tauInds_maxwell);
                                        if ~isempty(oldTau)
                                            ub_maxwell_loop(oldTau) = ub_loop_maxwell(oldTau);
                                            lb_maxwell_loop(oldTau) = lb_loop_maxwell(oldTau);
                                        end

                                        % Allow compliances to drop without restriction
                                        [oldMods,~] = ismember(find(~newInds_maxwell),modulusInds_maxwell);

                                        tauStep = floor(sqrt(loopLim));
                                        modStep = floor((loopLim)^(1./(length(ub_loop_maxwell(newMod)))))/(length(newTau)*tauStep);
                                        [betaGridE,betaGridTau] = meshgrid(logspace(log10(lb_loop_maxwell(newMod)),log10(ub_loop_maxwell(newMod)),modStep),...
                                            logspace(log10(lb_loop_maxwell(newTau)),log10(ub_loop_maxwell(newTau)),tauStep));

                                        if j <= tauStep*modStep
                                            beta0_maxwell_loop(newInds_maxwell) = [betaGridE(j);betaGridTau(j)];
                                        else
                                            beta0_maxwell_loop(newInds_maxwell) = [betaGridE(end);betaGridTau(end)];
                                        end

                                        beta0_maxwell_loop(oldMods) = beta0_maxwell_loop(oldMods)./oldModReductionFactor;

                                        % Enforce the original upper and lower
                                        % stiffness bounds
                                        if any(lb_maxwell_loop(oldMods) < lb_loop_maxwell(oldMods))
                                            lb_maxwell_loop(lb_maxwell_loop < lb_loop_maxwell) = lb_loop_maxwell(lb_maxwell_loop < lb_loop_maxwell);
                                        end
                                        if any(ub_maxwell_loop(oldMods) > ub_loop_maxwell(oldMods))
                                            ub_maxwell_loop(ub_maxwell_loop > ub_loop_maxwell) = ub_loop_maxwell(ub_maxwell_loop > ub_loop_maxwell);
                                        end

                                        if any(beta0_maxwell_loop > ub_maxwell_loop) || any(beta0_maxwell_loop < lb_maxwell_loop)
                                           badInds = logical((beta0_maxwell_loop > ub_maxwell_loop) + (beta0_maxwell_loop < lb_maxwell_loop));
                                           beta0_maxwell_loop(badInds) = (ub_maxwell_loop(badInds)-lb_maxwell_loop(badInds)) .* rand( size(beta0_maxwell_loop(badInds)) ) + lb_maxwell_loop(badInds);
                                        end

                                    else

                                        beta0_temp_loop = beta0_maxwell;
                                        if strcmp(elasticSetting_maxwell,'y')
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_maxwell(1)-lb_randLimit_maxwell(1))*rand()+lb_randLimit_maxwell(1)));
                                                beta0_temp_loop(2:2:end-1) = 10.^(((ub_randLimit_maxwell(2:2:end-1)-lb_randLimit_maxwell(2:2:end-1)).*rand(size(ub_randLimit_maxwell(2:2:end-1)))+lb_randLimit_maxwell(2:2:end-1)));
                                                beta0_temp_loop(3:2:end-1) = ((ub_randLimit_maxwell(3:2:end-1)-lb_randLimit_maxwell(3:2:end-1)).*rand(size(ub_randLimit_maxwell(3:2:end-1))) + lb_randLimit_maxwell(3:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_maxwell(1)-lb_randLimit_maxwell(1))*rand()+lb_randLimit_maxwell(1)));
                                                beta0_temp_loop(2:2:end) = 10.^(((ub_randLimit_maxwell(2:2:end)-lb_randLimit_maxwell(2:2:end)).*rand(size(ub_randLimit_maxwell(2:2:end)))+lb_randLimit_maxwell(2:2:end)));
                                                beta0_temp_loop(3:2:end) = ((ub_randLimit_maxwell(3:2:end)-lb_randLimit_maxwell(3:2:end)).*rand(size(ub_randLimit_maxwell(3:2:end))) + lb_randLimit_maxwell(3:2:end));
                                            end
                                        else
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1:2:end-1) = 10.^(((ub_randLimit_maxwell(1:2:end-1)-lb_randLimit_maxwell(1:2:end-1)).*rand(size(ub_randLimit_maxwell(1:2:end-1)))+lb_randLimit_maxwell(1:2:end-1)));
                                                beta0_temp_loop(2:2:end-1) = ((ub_randLimit_maxwell(2:2:end-1)-lb_randLimit_maxwell(2:2:end-1)).*rand(size(ub_randLimit_maxwell(2:2:end-1))) + lb_randLimit_maxwell(2:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1:2:end) = 10.^(((ub_randLimit_maxwell(1:2:end)-lb_randLimit_maxwell(1:2:end)).*rand(size(ub_randLimit_maxwell(1:2:end)))+lb_randLimit_maxwell(1:2:end)));
                                                beta0_temp_loop(2:2:end) = ((ub_randLimit_maxwell(2:2:end)-lb_randLimit_maxwell(2:2:end)).*rand(size(ub_randLimit_maxwell(2:2:end))) + lb_randLimit_maxwell(2:2:end));
                                            end
                                        end

                                        if strcmp(elasticSetting_maxwell,'y') && fitElasticFirst
                                            beta0_temp_loop(1) = 1/bestElasticTerm;
                                            ub_maxwell_loop = ub_loop_maxwell;
                                            lb_maxwell_loop = lb_loop_maxwell;
                                        else
                                            ub_maxwell_loop = ub_loop_maxwell;
                                            lb_maxwell_loop = lb_loop_maxwell;
                                        end

                                        beta0_maxwell_loop = beta0_temp_loop;

                                    end

                                    if any(newInds_schapery)
                                        beta0_temp_loop = beta0_schapery;
                                        if strcmp(elasticSetting_maxwell,'y')
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_schapery(1)-lb_randLimit_schapery(1))*rand()+lb_randLimit_schapery(1)));
                                                beta0_temp_loop(2:2:end-1) = 10.^(((ub_randLimit_schapery(2:2:end-1)-lb_randLimit_schapery(2:2:end-1)).*rand(size(ub_randLimit_schapery(2:2:end-1)))+lb_randLimit_schapery(2:2:end-1)));
                                                beta0_temp_loop(3:2:end-1) = ((ub_randLimit_schapery(3:2:end-1)-lb_randLimit_schapery(3:2:end-1)).*rand(size(ub_randLimit_schapery(3:2:end-1))) + lb_randLimit_schapery(3:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_schapery(1)-lb_randLimit_schapery(1))*rand()+lb_randLimit_schapery(1)));
                                                beta0_temp_loop(2:2:end) = 10.^(((ub_randLimit_schapery(2:2:end)-lb_randLimit_schapery(2:2:end)).*rand(size(ub_randLimit_schapery(2:2:end)))+lb_randLimit_schapery(2:2:end)));
                                                beta0_temp_loop(3:2:end) = ((ub_randLimit_schapery(3:2:end)-lb_randLimit_schapery(3:2:end)).*rand(size(ub_randLimit_schapery(3:2:end))) + lb_randLimit_schapery(3:2:end));
                                            end
                                        else
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1:2:end-1) = 10.^(((ub_randLimit_schapery(1:2:end-1)-lb_randLimit_schapery(1:2:end-1)).*rand(size(ub_randLimit_schapery(1:2:end-1)))+lb_randLimit_schapery(1:2:end-1)));
                                                beta0_temp_loop(2:2:end-1) = ((ub_randLimit_schapery(2:2:end-1)-lb_randLimit_schapery(2:2:end-1)).*rand(size(ub_randLimit_schapery(2:2:end-1))) + lb_randLimit_schapery(2:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1:2:end) = 10.^(((ub_randLimit_schapery(1:2:end)-lb_randLimit_schapery(1:2:end)).*rand(size(ub_randLimit_schapery(1:2:end)))+lb_randLimit_schapery(1:2:end)));
                                                beta0_temp_loop(2:2:end) = ((ub_randLimit_schapery(2:2:end)-lb_randLimit_schapery(2:2:end)).*rand(size(ub_randLimit_schapery(2:2:end))) + lb_randLimit_schapery(2:2:end));
                                            end
                                        end

                                        beta0_schapery_loop(newInds_schapery) = beta0_temp_loop(newInds_schapery);
                                        ub_schapery_loop = ub_loop_maxwell;
                                        lb_schapery_loop = lb_loop_maxwell;

                                        [oldMods,~] = ismember(find(~newInds_schapery),modulusInds_maxwell);
                                        beta0_schapery_loop(oldMods) = beta0_schapery_loop(oldMods)./oldModReductionFactor;

                                        % Enforce the original upper and lower
                                        % stiffness bounds
                                        if any(lb_schapery_loop(oldMods) < lb_loop_maxwell(oldMods))
                                            lb_schapery_loop(lb_schapery_loop < lb_loop_maxwell) = lb_loop_maxwell(lb_schapery_loop < lb_loop_maxwell);
                                        end
                                        if any(ub_schapery_loop(oldMods) > ub_loop_maxwell(oldMods))
                                            ub_schapery_loop(ub_schapery_loop > ub_loop_maxwell) = ub_loop_maxwell(ub_schapery_loop > ub_loop_maxwell);
                                        end

                                        % Enforce the original time
                                        % constant bounds
                                        [oldTau,~] = ismember(find(~newInds_schapery),tauInds_maxwell);
                                        if ~isempty(oldTau)
                                            ub_schapery_loop(oldTau) = ub_loop_maxwell(oldTau);
                                            lb_schapery_loop(oldTau) = lb_loop_maxwell(oldTau);
                                        end

                                        if any(beta0_schapery_loop > ub_maxwell_loop) || any(beta0_schapery_loop < lb_maxwell_loop)
                                           badInds = logical((beta0_schapery_loop > ub_maxwell_loop) + (beta0_schapery_loop < lb_maxwell_loop));
                                           beta0_schapery_loop(badInds) = (ub_maxwell_loop(badInds)-lb_maxwell_loop(badInds)) .* rand( size(beta0_schapery_loop(badInds)) ) + lb_maxwell_loop(badInds);
                                        end

                                    elseif any(newInds_schapery) && enforceGridGuess

                                        ub_schapery_loop = ub_loop_maxwell;
                                        lb_schapery_loop = lb_loop_maxwell;

                                        % Enforce the original time
                                        % constant bounds
                                        [oldTau,~] = ismember(find(~newInds_schapery),tauInds_maxwell);
                                        if ~isempty(oldTau)
                                            ub_schapery_loop(oldTau) = ub_loop_maxwell(oldTau);
                                            lb_schapery_loop(oldTau) = lb_loop_maxwell(oldTau);
                                        end

                                        % Allow compliances to drop without restriction
                                        [oldMods,~] = ismember(find(~newInds_schapery),modulusInds_maxwell);

                                        tauStep = floor(sqrt(loopLim));
                                        modStep = floor((loopLim)^(1./(length(ub_loop_maxwell(newMod)))))/(length(newTau)*tauStep);
                                        [betaGridE,betaGridTau] = meshgrid(logspace(log10(lb_loop_maxwell(newMod)),log10(ub_loop_maxwell(newMod)),modStep),...
                                            logspace(log10(lb_loop_maxwell(newTau)),log10(ub_loop_maxwell(newTau)),tauStep));

                                        if j <= tauStep*modStep
                                            beta0_schapery_loop(newInds_maxwell) = [betaGridE(j);betaGridTau(j)];
                                        else
                                            beta0_schapery_loop(newInds_maxwell) = [betaGridE(end);betaGridTau(end)];
                                        end

                                        beta0_schapery_loop(oldMods) = beta0_schapery_loop(oldMods)./oldModReductionFactor;

                                        % Enforce the original upper and lower
                                        % stiffness bounds
                                        if any(lb_schapery_loop(oldMods) < lb_loop_maxwell(oldMods))
                                            lb_schapery_loop(lb_schapery_loop < lb_loop_maxwell) = lb_loop_maxwell(lb_schapery_loop < lb_loop_maxwell);
                                        end
                                        if any(ub_schapery_loop(oldMods) > ub_loop_maxwell(oldMods))
                                            ub_schapery_loop(ub_schapery_loop > ub_loop_maxwell) = ub_loop_maxwell(ub_schapery_loop > ub_loop_maxwell);
                                        end

                                        if any(beta0_schapery_loop > ub_maxwell_loop) || any(beta0_schapery_loop < lb_maxwell_loop)
                                           badInds = logical((beta0_schapery_loop > ub_maxwell_loop) + (beta0_schapery_loop < lb_maxwell_loop));
                                           beta0_schapery_loop(badInds) = (ub_maxwell_loop(badInds)-lb_maxwell_loop(badInds)) .* rand( size(beta0_schapery_loop(badInds)) ) + lb_maxwell_loop(badInds);
                                        end

                                    else

                                        beta0_temp_loop = beta0_schapery;
                                        if strcmp(elasticSetting_maxwell,'y')
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_schapery(1)-lb_randLimit_schapery(1))*rand()+lb_randLimit_schapery(1)));
                                                beta0_temp_loop(2:2:end-1) = 10.^(((ub_randLimit_schapery(2:2:end-1)-lb_randLimit_schapery(2:2:end-1)).*rand(size(ub_randLimit_schapery(2:2:end-1)))+lb_randLimit_schapery(2:2:end-1)));
                                                beta0_temp_loop(3:2:end-1) = ((ub_randLimit_schapery(3:2:end-1)-lb_randLimit_schapery(3:2:end-1)).*rand(size(ub_randLimit_schapery(3:2:end-1))) + lb_randLimit_schapery(3:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1) = 10^(((ub_randLimit_schapery(1)-lb_randLimit_schapery(1))*rand()+lb_randLimit_schapery(1)));
                                                beta0_temp_loop(2:2:end) = 10.^(((ub_randLimit_schapery(2:2:end)-lb_randLimit_schapery(2:2:end)).*rand(size(ub_randLimit_schapery(2:2:end)))+lb_randLimit_schapery(2:2:end)));
                                                beta0_temp_loop(3:2:end) = ((ub_randLimit_schapery(3:2:end)-lb_randLimit_schapery(3:2:end)).*rand(size(ub_randLimit_schapery(3:2:end))) + lb_randLimit_schapery(3:2:end));
                                            end
                                        else
                                            if strcmp(fluidSetting_maxwell,'y')
                                                beta0_temp_loop(1:2:end-1) = 10.^(((ub_randLimit_schapery(1:2:end-1)-lb_randLimit_schapery(1:2:end-1)).*rand(size(ub_randLimit_schapery(1:2:end-1)))+lb_randLimit_schapery(1:2:end-1)));
                                                beta0_temp_loop(2:2:end-1) = ((ub_randLimit_schapery(2:2:end-1)-lb_randLimit_schapery(2:2:end-1)).*rand(size(ub_randLimit_schapery(2:2:end-1))) + lb_randLimit_schapery(2:2:end-1));
                                                beta0_temp_loop(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                                            else
                                                beta0_temp_loop(1:2:end) = 10.^(((ub_randLimit_schapery(1:2:end)-lb_randLimit_schapery(1:2:end)).*rand(size(ub_randLimit_schapery(1:2:end)))+lb_randLimit_schapery(1:2:end)));
                                                beta0_temp_loop(2:2:end) = ((ub_randLimit_schapery(2:2:end)-lb_randLimit_schapery(2:2:end)).*rand(size(ub_randLimit_schapery(2:2:end))) + lb_randLimit_schapery(2:2:end));
                                            end
                                        end

                                        if strcmp(elasticSetting_maxwell,'y') && fitElasticFirst
                                            beta0_temp_loop(1) = 1/bestElasticTerm/oldModReductionFactor;
                                            ub_schapery_loop = ub_loop_maxwell;
                                            lb_schapery_loop = lb_loop_maxwell;
                                        else
                                            ub_schapery_loop = ub_loop_maxwell;
                                            lb_schapery_loop = lb_loop_maxwell;
                                        end

                                        beta0_schapery_loop = beta0_temp_loop;

                                    end

                                    if plotCost
                                        legendEntries = {};
                                        for jj = 1:length(beta0_conv_loop)
                                            if ismember(jj,tauInds)
                                                if strcmp(fluidSetting,'y') && jj == length(beta0_conv_loop)
                                                    legendEntries = horzcat(legendEntries,sprintf('$\\phi_f$'));
                                                else
                                                    legendEntries = horzcat(legendEntries,sprintf('$\\tau_%d$',find(ismember(jj,tauInds))));
                                                end
                                            elseif ismember(jj,modulusInds)
                                                if strcmp(elasticSetting,'y') && jj == 1
                                                    legendEntries = horzcat(legendEntries,sprintf('$J_g$'));
                                                else
                                                    legendEntries = horzcat(legendEntries,sprintf('$J_%d$',find(ismember(jj,modulusInds))));
                                                end
                                            end
                                        end
                                    end

                                    % Gradient Descent
                                    if timeOperation
                                        Par.tic;
                                    end

                                    switch lower(optimizationMethod)
                                        case 'gen-alg'
                                            % Perform Genetic Algorithm Optimization
                                            opts = optimoptions('ga',...
                                                'FunctionTolerance',scaling);
                                            [newParams_voigt,finalRes,eflag,simoutput] = ga(@(x) sum((F_conv_loop(x,x_fit_loop)-y_fit_GD_loop).^2),length(beta0_conv_loop),[],[],[],[],lb_conv_loop,ub_conv_loop);

                                        case 'pattern'
                                            % Perform Pattern Search Optimization
                                            opts = optimoptions('patternsearch',...
                                                'FunctionTolerance',scaling);
                                            [newParams_voigt,finalRes,eflag,simoutput] = patternsearch(@(x) sum((F_conv_loop(x,x_fit_loop)-y_fit_GD_loop).^2),beta0_conv_loop,[],[],[],[],lb_conv_loop,ub_conv_loop);

                                        case 'particle-swarm'


                                        case 'surrogate'


                                    end

                                    if plotFit
                                        figure(fitPlot)
                                        scatter(x_fit_loop(:,1),y_fit_GD_loop,'bx')
                                        hold on
                                        plot(x_fit_loop(:,1),F_conv_loop(newParams_voigt,x_fit_loop),'r-')
                                        title('Final Optimized Model [Voigt]')
                                        xlabel('Time [s]')
                                        ylabel('Action Integral Value')
                                        set(gca,'YScale','log')
                                        set(gca,'XScale','log')
                                        legend('Data','Initial','Final')
                                        hold off
                                    end

                                    beta_dist(:,j) = newParams_voigt;
                                    beta0_dist(:,j) = zeros(size(beta0_conv));

                                    % Calculate Error Against Data Streams
                                    if forceResidual
                                        yModelIndentation = (1./alphaModel.*F_conv_loop(newParams_voigt,inputs_loop)).^(2/3);
                                        resnorm_dist(j) = sum(((yModelIndentation-ydata_voigt_loop).^2)./movvar(ydata_voigt_loop,3))./(length(ydata_voigt_loop)-length(newParams_voigt)); % Standard Error of Regression (S)
                                        voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',log_scale(yModelIndentation,x_fit_loop(:,1),dt,st)'];
                                    else
                                        resnorm_dist(j) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)
                                        voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',F_conv_wrapper_loop(newParams_voigt,x_fit_loop)'];
                                    end

                                    if timeOperation
                                        convTime(j) = Par.toc;
                                    end

                                    if plotCost
                                        figure(costPlot);
                                        clf;
                                    end

                                    if plotCost
                                        legendEntries = {};
                                        for jj = 1:length(beta0_maxwell_loop)
                                            if ismember(jj,tauInds_maxwell)
                                                if strcmp(fluidSetting_maxwell,'y') && jj == length(beta0_maxwell_loop)
                                                    legendEntries = horzcat(legendEntries,sprintf('$\\phi_f$'));
                                                else
                                                    legendEntries = horzcat(legendEntries,sprintf('$\\tau_%d$',find(ismember(jj,tauInds_maxwell))));
                                                end
                                            elseif ismember(jj,modulusInds_maxwell)
                                                if strcmp(elasticSetting_maxwell,'y') && jj == 1
                                                    legendEntries = horzcat(legendEntries,sprintf('$E_e$'));
                                                else
                                                    legendEntries = horzcat(legendEntries,sprintf('$E_%d$',find(ismember(jj,modulusInds_maxwell))));
                                                end
                                            end
                                        end
                                    end

                                    % Gradient Descent
                                    if timeOperation
                                        Par.tic;
                                    end

                                    switch lower(optimizationMethod)
                                        case 'gen-alg'
                                            % Perform Genetic Algorithm Optimization
                                            opts = optimoptions('ga',...
                                                'FunctionTolerance',scaling);
                                            [newParams_maxwell,finalRes,eflag,simoutput] = ga(@(x) sum((h_conv_schapery_area_loop(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2),length(beta0_maxwell_loop),[],[],[],[],lb_maxwell_loop,ub_maxwell_loop);

                                        case 'pattern'
                                            % Perform Pattern Search Optimization
                                            opts = optimoptions('patternsearch',...
                                                'FunctionTolerance',scaling);
                                            [newParams_maxwell,finalRes,eflag,simoutput] = patternsearch(@(x) sum((h_conv_schapery_area_loop(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2),beta0_maxwell_loop,[],[],[],[],lb_maxwell_loop,ub_maxwell_loop);

                                        case 'particle-swarm'


                                        case 'surrogate'


                                    end

                                    if plotFit
                                        figure(fitPlot)
                                        scatter(maxwellInputs_loop(:,1),y_fit_maxwell_GD_loop,'bx')
                                        hold on
                                        plot(maxwellInputs_loop(:,1),h_conv_maxwell_loop(newParams_maxwell,maxwellInputs_loop),'r-')
                                        title('Final Optimized Model [Maxwell]')
                                        xlabel('Time [s]')
                                        ylabel('Action Integral Value')
                                        set(gca,'YScale','log')
                                        set(gca,'XScale','log')
                                        legend('Data','Final')
                                        hold off
                                    end

                                    beta_dist_maxwell(:,j) = newParams_maxwell;
                                    beta0_dist_maxwell(:,j) = beta0_maxwell;

                                    % Calculate Error Against Data Streams
                                    if forceResidual
                                        if strcmp(elasticSetting_maxwell,'y')
                                            yModelForce = alphaModel.*(newParams_maxwell(1).*(inputsForce(:,2))+subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                        else
                                            yModelForce = alphaModel.*(subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                        end
                                        resnorm_dist_maxwell(j) = sum(((yModelForce-ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(newParams_maxwell)); % Standard Error of Regression (S)
                                        maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',log_scale(yModelForce,maxwellInputs_loop(:,1),dt,st)'];
                                    else
                                        resnorm_dist_maxwell(j) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)
                                        maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',h_conv_maxwell_wrapper_loop(newParams_maxwell,maxwellInputs_loop)'];
                                    end

                                    if timeOperation
                                        maxwellTime(j) = Par.toc;
                                    end

                                    if plotCost
                                        figure(costPlot);
                                        clf;
                                    end

                                    if plotCost
                                        legendEntries = {};
                                        for jj = 1:length(beta0_maxwell_loop)
                                            if ismember(jj,tauInds_maxwell)
                                                if strcmp(fluidSetting_maxwell,'y') && jj == length(beta0_maxwell_loop)
                                                    legendEntries = horzcat(legendEntries,sprintf('$\\phi_f$'));
                                                else
                                                    legendEntries = horzcat(legendEntries,sprintf('$\\tau_%d$',find(ismember(jj,tauInds_maxwell))));
                                                end
                                            elseif ismember(jj,modulusInds_maxwell)
                                                if strcmp(elasticSetting_maxwell,'y') && jj == 1
                                                    legendEntries = horzcat(legendEntries,sprintf('$E_e$'));
                                                else
                                                    legendEntries = horzcat(legendEntries,sprintf('$E_%d$',find(ismember(jj,modulusInds_maxwell))));
                                                end
                                            end
                                        end
                                    end

                                    % Gradient Descent
                                    if timeOperation
                                        Par.tic;
                                    end

                                    switch lower(optimizationMethod)
                                        case 'gen-alg'
                                            % Perform Genetic Algorithm Optimization
                                            opts = optimoptions('ga',...
                                                'FunctionTolerance',scaling);
                                            [newParams_schapery,finalRes,eflag,simoutput] = ga(@(x) sum((h_conv_schapery_area_loop(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2),length(beta0_schapery_loop),[],[],[],[],lb_schapery_loop,ub_schapery_loop);

                                        case 'pattern'
                                            % Perform Pattern Search Optimization
                                            opts = optimoptions('patternsearch',...
                                                'FunctionTolerance',scaling);
                                            [newParams_schapery,finalRes,eflag,simoutput] = patternsearch(@(x) sum((h_conv_schapery_area_loop(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2),beta0_schapery_loop,[],[],[],[],lb_schapery_loop,ub_schapery_loop);

                                        case 'particle-swarm'


                                        case 'surrogate'


                                    end

                                    if plotFit
                                        figure(fitPlot)
                                        scatter(x_fit_area_loop(:,1),y_fit_schapery_GD_loop,'bx')
                                        hold on
                                        plot(x_fit_area_loop(:,1),h_conv_schapery_area_loop(newParams_schapery,x_fit_area_loop),'r-')
                                        title('Final Optimized Model [Schapery-Maxwell]')
                                        xlabel('Time [s]')
                                        ylabel('Force [N]')
                                        set(gca,'YScale','log')
                                        set(gca,'XScale','log')
                                        legend('Data','Final')
                                        hold off
                                    end

                                    beta_dist_schapery_area(:,j) = newParams_schapery;
                                    beta0_dist_schapery_area(:,j) = zeros(size(beta0_schapery));

                                    % Calculate Error Against Data Streams
                                    if forceResidual
                                        yModelForce = h_conv_schapery_area_loop(newParams_schapery,inputsForce_schapery);
                                        resnorm_dist_schapery_area(j) = sum(((yModelForce-ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(newParams_schapery)); % Standard Error of Regression (S)
                                        schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',log_scale(yModelForce,x_fit_area_loop(:,1),dt,st)'];
                                    else
                                        resnorm_dist_schapery_area(j) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)
                                        schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',h_conv_schapery_area_wrapper_loop(newParams_schapery,x_fit_area_loop)'];
                                    end

                                    if timeOperation
                                        schaperyTime(j) = Par.toc;
                                    end

                                    if timeOperation
                                        % Track memory
                                        %ramused(j) = userMem.MemUsedMATLAB;
                                        ramused(j) = 0;
                                        %ramav(j) = systemMem.PhysicalMemory.Available;
                                        ramav(j) = 0;
                                        ramtime(j) = now;
                                    end

                                    % Placeholder for the linear model since we
                                    % don't care about it.
                                    beta_dist_linear(:,j) = ones(size(beta0_conv));
                                    beta0_dist_linear(:,j) = ones(size(beta0_conv));
                                    resnorm_dist_linear(j) = 1;

                                end % End Parfor

                            case 'nelder-mead'

                                ydata_loop = F_r_loop;
                                ydata_voigt_loop = h_r_loop;
                                inputs_loop = [x_fit_loop(:,1),...
                                        F_r_loop];
                                inputsForce = [x_fit_loop(:,1),...
                                    h_r_loop.^1.5];
                                inputsForce_schapery = [x_fit_loop(x_fit_loop(:,1)>0,1),...
                                    sqrt(2*r_tip*h_r_loop(x_fit_loop(:,1)>0)...
                                    - (h_r_loop(x_fit_loop(:,1)>0).^2))./r_tip];
                                timeVec_loop = x_fit_loop(:,1);
                                alphaModel = (8*sqrt(r_tip))/(3*(1-nu_sample));
                                y_fit_GD_loop = h_norm_loop;
                                y_fit_maxwell_GD_loop = F_norm_loop;
                                y_fit_schapery_GD_loop = F_schapery_area_loop;

                                useAnnealing = 0;
                                tossOutOfBounds = 0;

                                parfor j = 1:loopLim
                                    warning('off');
                                    % Initialize
                                    newInds_linear = [];
                                    newInds_conv = [];
                                    newInds_maxwell = [];
                                    newInds_schapery = [];
                                    beta0_linear_loop = [];
                                    beta0_conv_loop = [];
                                    beta0_maxwell_loop = [];
                                    beta0_schapery_loop = [];
                                    ub_randLimit = [];
                                    lb_randLimit = [];
                                    ub_randLimit_maxwell = [];
                                    lb_randLimit_maxwell = [];
                                    ub_randLimit_schapery = [];
                                    lb_randLimit_schapery = [];
                                    beta0_temp_loop = [];
                                    limScale = [];
                                    tauInds = [];
                                    oldTau = [];
                                    modulusInds = [];
                                    oldModuli = [];
                                    newParams_temp = [];
                                    F_conv_loop_weight = [];
                                    h_conv_maxwell_loop_weight = [];
                                    h_conv_schapery_area_loop_weight = [];
                                    beta0_linear_input = [];
                                    beta0_linear_input = beta0_linear;
                                    beta0_conv_input = [];
                                    beta0_conv_input = beta0_conv;
                                    beta0_maxwell_input = [];
                                    beta0_maxwell_input = beta0_maxwell;
                                    beta0_schapery_input = [];
                                    beta0_schapery_input = beta0_schapery;

                                    % Check for bad indices being fed forward
                                    if ( any((ub_loop(~isnan(beta0_linear)) - beta0_linear(~isnan(beta0_linear))) < 0) ) || ( any((beta0_linear(~isnan(beta0_linear)) - lb_loop(~isnan(beta0_linear))) < 0) )
                                        beta0_linear_input(or((ub_loop(~isnan(beta0_linear)) - beta0_linear(~isnan(beta0_linear))) < 0,(beta0_linear(~isnan(beta0_linear)) - lb_loop(~isnan(beta0_linear))) < 0)) = NaN;
                                    end
                                    if ( any((ub_loop(~isnan(beta0_conv)) - beta0_conv(~isnan(beta0_conv))) < 0) ) || ( any((beta0_conv(~isnan(beta0_conv)) - lb_loop(~isnan(beta0_conv))) < 0) )
                                        beta0_conv_input(or((ub_loop(~isnan(beta0_conv)) - beta0_conv(~isnan(beta0_conv))) < 0,(beta0_conv(~isnan(beta0_conv)) - lb_loop(~isnan(beta0_conv))) < 0)) = NaN;
                                    end
                                    if ( any((ub_loop_maxwell(~isnan(beta0_maxwell)) - beta0_maxwell(~isnan(beta0_maxwell))) < 0) ) || ( any((beta0_maxwell(~isnan(beta0_maxwell)) - lb_loop_maxwell(~isnan(beta0_maxwell))) < 0) )
                                        beta0_maxwell_input(or((ub_loop_maxwell(~isnan(beta0_maxwell)) - beta0_maxwell(~isnan(beta0_maxwell))) < 0,(beta0_maxwell(~isnan(beta0_maxwell)) - lb_loop_maxwell(~isnan(beta0_maxwell))) < 0)) = NaN;
                                    end
                                    if ( any((ub_loop_maxwell(~isnan(beta0_schapery)) - beta0_schapery(~isnan(beta0_schapery))) < 0) ) || ( any((beta0_schapery(~isnan(beta0_schapery)) - lb_loop_maxwell(~isnan(beta0_schapery))) < 0) )
                                        beta0_schapery_input(or((ub_loop_maxwell(~isnan(beta0_schapery)) - beta0_schapery(~isnan(beta0_schapery))) < 0,(beta0_schapery(~isnan(beta0_schapery)) - lb_loop_maxwell(~isnan(beta0_schapery))) < 0)) = NaN;
                                    end

                                    % Find the new indices
                                    newInds_linear = isnan(beta0_linear_input);
                                    newInds_conv = isnan(beta0_conv_input);
                                    newInds_maxwell = isnan(beta0_maxwell_input);
                                    newInds_schapery = isnan(beta0_schapery_input);

                                    % Control the bounds AND the initial parameters
                                    ub_randLimit = ub_randLimit_multigrid;
                                    lb_randLimit = lb_randLimit_multigrid;
                                    ub_randLimit_maxwell = ub_randLimit_maxwell_multigrid;
                                    lb_randLimit_maxwell = lb_randLimit_maxwell_multigrid;
                                    ub_randLimit_schapery = ub_randLimit_schapery_multigrid;
                                    lb_randLimit_schapery = lb_randLimit_schapery_multigrid;
                                    ub_fluidityRandLimit_init = 1e-20;
                                    limScale = getfield(logspace(5,0,n_terms),{i})*relaxationFactor;

                                    % Include or remove elastic term from bound
                                    % loosening that happens below. Use 1 for
                                    % inclusion, empty brackets [] for ignoring.
                                    elasticInd = 1;

        %                             oldModReductionFactor = limScale*0.5;
                                    oldModReductionFactor = 1;
        %                             oldModReductionFactor = 1/relaxationFactor;

                                    options = optimset('Display','none',...
                                        'PlotFcns',[],...
                                        'MaxFunEvals',n_maxIterations,...
                                        'MaxIter',n_maxIterations,...
                                        'TolFun',scaling,...
                                        'TolX',0);

                                    nelderopts = struct(...
                                        'CoolSched',@(T) (.8*T),...
                                        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
                                        'InitTemp',1,...
                                        'MaxConsRej',1000,...
                                        'MaxSuccess',20,...
                                        'MaxTries',n_maxIterations,...
                                        'StopTemp',scaling,...
                                        'StopVal',-Inf,...
                                        'Verbosity',0);

                                    % Create random starting point
                                    [beta0_conv_loop,tauInds,modulusInds] = makeRandomParams(beta0_conv_input,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,newInds_conv);

                                    weightArray = ones(size(beta0_conv_loop));
                                    if scaleParams
                                        numTau = length(tauInds);
                                        if strcmp(fluidSetting,'y')
                                           numTau = numTau-1;
                                        end
                                        if i == 1
                                            switch lower(paramScale)
                                                case 'logscale'
                                                    F_conv_loop_weight = logscaleParams(F_conv_loop);

                                                otherwise
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(modulusInds) = 10.^((max(ub_randLimit(modulusInds))-min(lb_randLimit(modulusInds)))/2);
                                                    F_conv_loop_weight = weightParams(F_conv_loop,weightArray);
                                            end
                                        else
                                            switch lower(paramScale)
                                                case 'equal'
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_conv_input(intersect(find(~newInds_conv),modulusInds))))));
                                                    F_conv_loop_weight = weightParams(F_conv_loop,weightArray);

                                                case 'individual'
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(intersect(find(~newInds_conv),modulusInds)) = 10.^(floor(log10(beta0_conv_input(intersect(find(~newInds_conv),modulusInds)))));
                                                    weightArray(intersect(find(newInds_conv),modulusInds)) = mean(10.^(floor(log10(beta0_conv_input(intersect(find(~newInds_conv),modulusInds))))));
                                                    F_conv_loop_weight = weightParams(F_conv_loop,weightArray);

                                                case 'logscale'
                                                    F_conv_loop_weight = logscaleParams(F_conv_loop);
                                            end
                                        end
                                    else
                                        F_conv_loop_weight = weightParams(F_conv_loop,weightArray);
                                    end

                                    if scaleParams
                                        switch lower(paramScale)
                                            case 'logscale'
                                                ub_conv_loop = log10(ub_loop);
                                                lb_conv_loop = log10(lb_loop);
                                                beta0_conv_loop = log10(beta0_conv_loop);
                                            otherwise
                                                ub_conv_loop = (ub_loop)./weightArray;
                                                lb_conv_loop = (lb_loop)./weightArray;
                                                beta0_conv_loop = beta0_conv_loop./weightArray;
                                        end
                                    else
                                        ub_conv_loop = (ub_loop);
                                        lb_conv_loop = (lb_loop);
                                    end

                                    % Gradient Descent
                                    if timeOperation
                                        Par.tic;
                                    end

                                    costObj = @(x) sum((F_conv_loop_weight(x,x_fit_loop)-y_fit_GD_loop).^2);

                                    if ~useAnnealing
                                        [newParams_voigt,finalRes,eflag] = fminsearch(@(c) boundObjective(costObj,c,ub_conv_loop,lb_conv_loop), beta0_conv_loop, options);
                                    else
                                        [newParams_voigt,finalRes] = annealOpt(@(c) boundObjective(costObj,c,ub_conv_loop,lb_conv_loop), beta0_conv_loop, nelderopts, options);
                                    end

                                    if ~strcmp(lower(paramScale),'logscale') && scaleParams
                                        newParams_voigt = newParams_voigt.*weightArray;
                                        beta0_conv_loop = beta0_conv_loop.*weightArray;
                                    elseif  strcmp(lower(paramScale),'logscale') && scaleParams
                                        newParams_voigt = 10.^(newParams_voigt);
                                        beta0_conv_loop = 10.^(beta0_conv_loop);
                                    end

                                    if tossOutOfBounds
                                        if ( any((ub_loop - newParams_voigt) < 0) ) || ( any((newParams_voigt - lb_loop) < 0) )
                                            newParams_voigt = NaN(size(newParams_voigt));
                                        end
                                    else
                                        if any(newParams_voigt < 0)
                                            newParams_voigt = NaN(size(newParams_voigt));
                                        end
                                    end

                                    beta_dist(:,j) = newParams_voigt;
                                    beta0_dist(:,j) = beta0_conv_loop;

                                    % Calculate Error Against Data Streams
                                    if forceResidual
                                        yModelIndentation = (1./alphaModel.*F_conv_loop(newParams_voigt,inputs_loop)).^(2/3);
                                        resnorm_dist(j) = sum(((yModelIndentation-ydata_voigt_loop).^2)./movvar(ydata_voigt_loop,3))./(length(ydata_voigt_loop)-length(newParams_voigt)); % Standard Error of Regression (S)
                                        voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',log_scale(yModelIndentation,x_fit_loop(:,1),dt,st)'];
                                    else
                                        resnorm_dist(j) = errorCalcVoigt(newParams_voigt,x_fit_loop); % Standard Error of Regression (S)
                                        voigtData = [log_scale(x_fit_loop(:,1),x_fit_loop(:,1),dt,st)',F_conv_wrapper_loop(newParams_voigt,x_fit_loop)'];
                                    end

                                    if timeOperation
                                        convTime(j) = Par.toc;
                                    end

                                    % Create random starting point
                                    [beta0_maxwell_loop,tauInds,modulusInds] = makeRandomParams(beta0_maxwell_input,ub_randLimit_maxwell,lb_randLimit_maxwell,elasticSetting_maxwell,fluidSetting_maxwell,newInds_maxwell);

                                    weightArray = ones(size(beta0_maxwell_loop));
                                    if scaleParams
                                        numTau = length(tauInds);
                                        if strcmp(fluidSetting_maxwell,'y')
                                           numTau = numTau-1;
                                        end
                                        if i == 1
                                            switch lower(paramScale)
                                                case 'logscale'
                                                    h_conv_maxwell_loop_weight = logscaleParams(h_conv_maxwell_loop);

                                                otherwise
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(modulusInds) = 10.^((max(ub_randLimit_maxwell(modulusInds))-min(lb_randLimit_maxwell(modulusInds)))/2);
                                                    h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);

                                            end
                                        else
                                            switch lower(paramScale)
                                                case 'equal'
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_maxwell_input(~newInds_maxwell)))));
                                                    h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);

                                                case 'individual'
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(intersect(find(~newInds_maxwell),modulusInds)) = 10.^(floor(log10(beta0_maxwell_input(intersect(find(~newInds_maxwell),modulusInds)))));
                                                    weightArray(intersect(find(newInds_maxwell),modulusInds)) = mean(10.^(floor(log10(beta0_maxwell_input(intersect(find(~newInds_maxwell),modulusInds))))));
                                                    h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);

                                                case 'logscale'
                                                    h_conv_maxwell_loop_weight = logscaleParams(h_conv_maxwell_loop);

                                            end
                                        end
                                    else
                                        h_conv_maxwell_loop_weight = weightParams(h_conv_maxwell_loop,weightArray);
                                    end

                                    if scaleParams
                                        switch lower(paramScale)
                                            case 'logscale'
                                                ub_maxwell_loop = log10(ub_loop_maxwell);
                                                lb_maxwell_loop = log10(lb_loop_maxwell);
                                                beta0_maxwell_loop = log10(beta0_maxwell_loop);
                                            otherwise
                                                ub_maxwell_loop = (ub_loop_maxwell)./weightArray;
                                                lb_maxwell_loop = (lb_loop_maxwell)./weightArray;
                                                beta0_maxwell_loop = beta0_maxwell_loop./weightArray;
                                        end
                                    else
                                        ub_maxwell_loop = (ub_loop_maxwell);
                                        lb_maxwell_loop = (lb_loop_maxwell);
                                    end

                                    % Gradient Descent
                                    if timeOperation
                                        Par.tic;
                                    end

                                    costObj = @(x) sum((h_conv_maxwell_loop_weight(x,maxwellInputs_loop)-y_fit_maxwell_GD_loop).^2);

                                    if ~useAnnealing
                                        [newParams_maxwell,finalRes,eflag] = fminsearch(@(c) boundObjective(costObj,c,ub_maxwell_loop,lb_maxwell_loop), beta0_maxwell_loop, options);
                                    else
                                        [newParams_maxwell,finalRes] = annealOpt(@(c) boundObjective(costObj,c,ub_maxwell_loop,lb_maxwell_loop), beta0_maxwell_loop, nelderopts, options);
                                    end

                                    if ~strcmp(lower(paramScale),'logscale') && scaleParams
                                        newParams_maxwell = newParams_maxwell.*weightArray;
                                        beta0_maxwell_loop = beta0_maxwell_loop.*weightArray;
                                    elseif strcmp(lower(paramScale),'logscale') && scaleParams
                                        newParams_maxwell = 10.^(newParams_maxwell);
                                        beta0_maxwell_loop = 10.^(beta0_maxwell_loop);
                                    end

                                    if tossOutOfBounds
                                        if ( any((ub_loop_maxwell - newParams_maxwell) < 0) ) || ( any((newParams_maxwell - lb_loop_maxwell) < 0) )
                                            newParams_maxwell = NaN(size(newParams_maxwell));
                                        end
                                    else
                                        if any(newParams_maxwell < 0)
                                            newParams_maxwell = NaN(size(newParams_maxwell));
                                        end
                                    end

                                    beta_dist_maxwell(:,j) = newParams_maxwell;
                                    beta0_dist_maxwell(:,j) = beta0_maxwell_loop;

                                    % Calculate Error Against Data Streams
                                    if forceResidual
                                        if strcmp(elasticSetting_maxwell,'y')
                                            yModelForce = alphaModel.*(newParams_maxwell(1).*(inputsForce(:,2))+subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                        else
                                            yModelForce = alphaModel.*(subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce), gradient(inputsForce(:,2),dt),'full'))*dt);
                                        end
                                        resnorm_dist_maxwell(j) = sum(((yModelForce-ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(newParams_maxwell)); % Standard Error of Regression (S)
                                        maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',log_scale(yModelForce,maxwellInputs_loop(:,1),dt,st)'];
                                    else
                                        resnorm_dist_maxwell(j) = errorCalcMaxwell(newParams_maxwell,maxwellInputs_loop); % Standard Error of Regression (S)
                                        maxwellData = [log_scale(maxwellInputs_loop(:,1),maxwellInputs_loop(:,1),dt,st)',h_conv_maxwell_wrapper_loop(newParams_maxwell,maxwellInputs_loop)'];
                                    end

                                    if timeOperation
                                        maxwellTime(j) = Par.toc;
                                    end

                                    % Create random starting point
                                    [beta0_schapery_loop,tauInds,modulusInds] = makeRandomParams(beta0_schapery_input,ub_randLimit_schapery,lb_randLimit_schapery,elasticSetting_maxwell,fluidSetting_maxwell,newInds_schapery);

                                    weightArray = ones(size(beta0_schapery_loop));
                                    if scaleParams
                                        numTau = length(tauInds);
                                        if strcmp(fluidSetting_maxwell,'y')
                                           numTau = numTau-1;
                                        end
                                        if i == 1
                                            switch lower(paramScale)
                                                case 'logscale'
                                                    h_conv_schapery_area_loop_weight = logscaleParams(h_conv_schapery_area_loop);

                                                otherwise
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(modulusInds) = 10.^((max(ub_randLimit_schapery(modulusInds))-min(lb_randLimit_schapery(modulusInds)))/2);
                                                    h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);

                                            end
                                        else
                                            switch lower(paramScale)
                                                case 'equal'
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_schapery_input(~newInds_schapery)))));
                                                    h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);

                                                case 'individual'
                                                    weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                    weightArray(intersect(find(~newInds_schapery),modulusInds)) = 10.^(floor(log10(beta0_schapery_input(intersect(find(~newInds_schapery),modulusInds)))));
                                                    weightArray(intersect(find(newInds_schapery),modulusInds)) = mean(10.^(floor(log10(beta0_schapery_input(intersect(find(~newInds_schapery),modulusInds))))));
                                                    h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);

                                                case 'logscale'
                                                    h_conv_schapery_area_loop_weight = logscaleParams(h_conv_schapery_area_loop);

                                            end
                                        end
                                    else
                                        h_conv_schapery_area_loop_weight = weightParams(h_conv_schapery_area_loop,weightArray);
                                    end

                                    if scaleParams
                                        switch lower(paramScale)
                                            case 'logscale'
                                                ub_schapery_loop = log10(ub_loop_maxwell);
                                                lb_schapery_loop = log10(lb_loop_maxwell);
                                                beta0_schapery_loop = log10(beta0_schapery_loop);
                                            otherwise
                                                ub_schapery_loop = (ub_loop_maxwell)./weightArray;
                                                lb_schapery_loop = (lb_loop_maxwell)./weightArray;
                                                beta0_schapery_loop = beta0_schapery_loop./weightArray;
                                        end
                                    else
                                        ub_schapery_loop = (ub_loop_maxwell);
                                        lb_schapery_loop = (lb_loop_maxwell);
                                    end

                                    % Gradient Descent
                                    if timeOperation
                                        Par.tic;
                                    end

                                    costObj = @(x) sum((h_conv_schapery_area_loop_weight(x,x_fit_area_loop)-y_fit_schapery_GD_loop).^2);

                                    if ~useAnnealing
                                        [newParams_schapery,finalRes,eflag] = fminsearch(@(c) boundObjective(costObj,c,ub_schapery_loop,lb_schapery_loop), beta0_schapery_loop, options);
                                    else
                                        [newParams_schapery,finalRes] = annealOpt(@(c) boundObjective(costObj,c,ub_schapery_loop,lb_schapery_loop), beta0_schapery_loop, nelderopts, options);
                                    end

                                    if ~strcmp(lower(paramScale),'logscale') && scaleParams
                                        newParams_schapery = newParams_schapery.*weightArray;
                                        beta0_schapery_loop = beta0_schapery_loop.*weightArray;
                                    elseif strcmp(lower(paramScale),'logscale') && scaleParams
                                        newParams_schapery = 10.^(newParams_schapery);
                                        beta0_schapery_loop = 10.^(beta0_schapery_loop);
                                    end

                                    if tossOutOfBounds
                                        if ( any((ub_loop_maxwell - newParams_schapery) < 0) ) || ( any((newParams_schapery - lb_loop_maxwell) < 0) )
                                            newParams_schapery = NaN(size(newParams_schapery));
                                        end
                                    else
                                        if any(newParams_schapery < 0)
                                            newParams_schapery = NaN(size(newParams_schapery));
                                        end
                                    end

                                    beta_dist_schapery_area(:,j) = newParams_schapery;
                                    beta0_dist_schapery_area(:,j) = beta0_schapery_loop;

                                    % Calculate Error Against Data Streams
                                    if forceResidual
                                        yModelForce = h_conv_schapery_area_loop(newParams_schapery,inputsForce_schapery);
                                        resnorm_dist_schapery_area(j) = sum(((yModelForce-ydata_loop).^2)./movvar(ydata_loop,3))./(length(ydata_loop)-length(newParams_schapery)); % Standard Error of Regression (S)
                                        schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',log_scale(yModelForce,x_fit_area_loop(:,1),dt,st)'];
                                    else
                                        resnorm_dist_schapery_area(j) = errorCalcSchapery(newParams_schapery,x_fit_area_loop); % Standard Error of Regression (S)
                                        schaperyData = [log_scale(x_fit_area_loop(:,1),x_fit_area_loop(:,1),dt,st)',h_conv_schapery_area_wrapper_loop(newParams_schapery,x_fit_area_loop)'];
                                    end

                                    if timeOperation
                                        schaperyTime(j) = Par.toc;
                                    end

                                    if timeOperation
                                        % Track memory
                                        %ramused(j) = userMem.MemUsedMATLAB;
                                        ramused(j) = 0;
                                        %ramav(j) = systemMem.PhysicalMemory.Available;
                                        ramav(j) = 0;
                                        ramtime(j) = now;
                                    end

                                    % Placeholder for the linear model since we
                                    % don't care about it.
                                    beta_dist_linear(:,j) = ones(size(beta0_conv_loop));
                                    beta0_dist_linear(:,j) = ones(size(beta0_conv_loop));
                                    resnorm_dist_linear(j) = 1;

                                end % End Parfor

                            otherwise
                                error('That optimization method is not supported for the iterative approach.')
                        end

                        if timeOperation
                            parallelFitEndTime(iVoigt,i_loop,i) = toc;
                            parallelFitParallelBytes(iVoigt,i_loop,i) = getfield(tocBytes(gcp),{1});

                            stop(convTime);
                            stop(maxwellTime);
                            stop(schaperyTime);

                            convTimeTotal = [convTimeTotal convTime];
                            maxwellTimeTotal = [maxwellTimeTotal maxwellTime];
                            schaperyTimeTotal = [schaperyTimeTotal schaperyTime];

                            fprintf('\nThe parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                                (parallelFitEndTime(iVoigt,i_loop,i)-parallelFitStartTime(iVoigt,i_loop,i))/60, parallelFitParallelBytes(iVoigt,i_loop,i)*1e-6);
                        end

                        resnorm_dist_linear_temp = resnorm_dist_linear;
                        resnorm_dist_linear_temp(resnorm_dist_linear_temp == 0) = NaN;
                        [~,idx_linear] = min(resnorm_dist_linear_temp,[],'omitnan');
                        linearParams_loop{i} = beta_dist_linear(:,idx_linear);

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

                        resnorm_dist_schapery_area_temp = resnorm_dist_schapery_area;
                        resnorm_dist_schapery_area_temp(resnorm_dist_schapery_area_temp == 0) = NaN;
                        resnorm_dist_schapery_area_temp( any(abs(log10(beta_dist_schapery_area(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                        resnorm_dist_schapery_area_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))) - log10(beta_dist_schapery_area(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                        [~,idx_schapery] = min(resnorm_dist_schapery_area_temp,[],'omitnan');
                        schaperyParams_loop{i} = beta_dist_schapery_area(:,idx_schapery);

                    end

                    if timeOperation
                        iterativeFitEndTime(iVoigt,i_loop) = toc;
                        parallelFitBreakdown{iVoigt,i_loop} = {convTimeTotal,maxwellTimeTotal,schaperyTimeTotal};
                        parallelFitMemoryMonitor{iVoigt,i_loop} = vertcat(ramused,ramav,ramtime);
                    end

                    dataStruct(indShift+i_loop).linearParams_loop = linearParams_loop;
                    dataStruct(indShift+i_loop).convParams_loop = convParams_loop;
                    dataStruct(indShift+i_loop).schaperyParams_loop = schaperyParams_loop;
                    dataStruct(indShift+i_loop).maxwellParams_loop = maxwellParams_loop;

                    n_ranks = 10; % Number of places to save the parameter sets of

    %                 if n_ranks < length(beta_dist(1,~isnan(beta_dist(1,:))))
    %                     n_ranks = length(beta_dist(1,~isnan(beta_dist(1,:))));
    %                 end

    %                 [~,idx_linear] = min(resnorm_dist_linear);
    %                 dataStruct(indShift+i_loop).linearParams = beta_dist_linear(:,idx_linear);
    % 
    %                 [~,idx] = min(resnorm_dist(resnorm_dist > 0));
    %                 dataStruct(indShift+i_loop).convParams = beta_dist(:,idx);
    %                 
    %                 [~,I] = sort(resnorm_dist(resnorm_dist > 0));
    %                 dataStruct(indShift+i_loop).convParams_alternatives = beta_dist(:,I(1:n_ranks));
    %                 
    %                 [~,idx_maxwell] = min(resnorm_dist_maxwell);
    %                 dataStruct(indShift+i_loop).maxwellParams = beta_dist_maxwell(:,idx_maxwell);
    % 
    %                 [~,I_maxwell] = sort(resnorm_dist_maxwell);
    %                 dataStruct(indShift+i_loop).maxwellParams_alternatives = beta_dist_maxwell(:,I_maxwell(1:n_ranks));
    % 
    %                 [~,idx_schapery] = min(resnorm_dist_schapery_area);
    %                 dataStruct(indShift+i_loop).schaperyParams = beta_dist_schapery_area(:,idx_schapery);
    %                 
    %                 [~,I_schapery] = sort(resnorm_dist_schapery_area);
    %                 dataStruct(indShift+i_loop).schaperyParams_alternatives = beta_dist_schapery_area(:,I_schapery(1:n_ranks)); 

                    resnorm_dist_linear_temp = resnorm_dist_linear;
                    resnorm_dist_linear_temp(resnorm_dist_linear_temp == 0) = NaN;
                    [~,idx_linear] = min(resnorm_dist_linear_temp,[],'omitnan');
                    dataStruct(indShift+i_loop).linearParams = beta_dist_linear(:,idx_linear);

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

                    resnorm_dist_schapery_area_temp = resnorm_dist_schapery_area;
                    resnorm_dist_schapery_area_temp(resnorm_dist_schapery_area_temp == 0) = NaN;
                    resnorm_dist_schapery_area_temp( any(abs(log10(beta_dist_schapery_area(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                    resnorm_dist_schapery_area_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))) - log10(beta_dist_schapery_area(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                    [~,idx_schapery] = min(resnorm_dist_schapery_area_temp,[],'omitnan');
                    dataStruct(indShift+i_loop).schaperyParams = beta_dist_schapery_area(:,idx_schapery);

                    [~,I_schapery] = sort(resnorm_dist_schapery_area_temp,'MissingPlacement','last');
                    dataStruct(indShift+i_loop).schaperyParams_alternatives = beta_dist_schapery_area(:,I_schapery(1:n_ranks));

    %                 if i == n_terms
    %                     disp('Paused on Last Results');
    %                     % Searching to see if the RIGHT parameters even exist
    %                     % in the set!
    %                 end

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
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,idx_maxwell),iii,beta_dist_maxwell(iii*2,idx_maxwell));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                    fprintf(fid,'\r\n\r\nSchapery Generalized Maxwell\r\n');
                    if elasticSetting_maxwell == 'y'
                        fprintf(fid,'E_e = %3.6g [Pa]\r\n',beta_dist_schapery_area(1,idx_schapery));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,idx_schapery),iii,beta_dist_schapery_area(iii*2+1,idx_schapery));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,idx_schapery),iii,beta_dist_schapery_area(iii*2,idx_schapery));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
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
                            fprintf(fid,'E_g = %3.6g [Pa]\r\n',beta_dist_maxwell(1,I_maxwell(alti)));

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2+1,I_maxwell(alti)));
                            end
                            if fluidSetting_maxwell == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        else
                            fprintf(fid,'No elastic response used.\r\n');

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2,I_maxwell(alti)));
                            end
                            if fluidSetting_maxwell == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        end

                        % Print Schapery Maxwell
                        fprintf(fid,'\r\nSchapery Generalized Maxwell\n');
                        if elasticSetting_maxwell == 'y'
                            fprintf(fid,'E_g = %3.6g [Pa]\r\n',beta_dist_schapery_area(1,I_schapery(alti)));

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2+1,I_schapery(alti)));
                            end
                            if fluidSetting_maxwell == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        else
                            fprintf(fid,'No elastic response used.\r\n');

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2,I_schapery(alti)));
                            end
                            if fluidSetting_maxwell == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        end

                    end

                    % Calculate Harmonic Quantities
                    if smoothData == 'n'
                        inputs = [dataStruct(indShift+i_loop).t_r, dataStruct(indShift+i_loop).F_r];
                    else
                        inputs = [dataStruct(indShift+i_loop).t_r_smooth, dataStruct(indShift+i_loop).F_r_smooth];
                    end

                    % Define the convolution fit dataset generation function
                    data_fit = cumtrapz(F_conv(dataStruct(indShift+i_loop).convParams,inputs)).*dt;
                    data_fit_log = log_scale(data_fit, timeVec, dt, st);

                    % Define the linear convolution fit dataset generation function
                    linear_fit = cumtrapz(F_conv(dataStruct(indShift+i_loop).linearParams,inputs)).*dt;
                    linear_fit_log = log_scale(linear_fit, timeVec, dt, st);

                    % Move on to plotting the loss angle
                    % First, generate a frequency array in log scale
                    de0 = 2.0.*pi.*(1./timeVec(end));
                    maxi = 2.0.*pi.*(1./mean(diff(timeVec)));
                    omega = log_tw(de0,maxi);

                    % Second, use the convolution J parameters
                    % Calculate J_storage
                    Jstorage = J_storage_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting);

                    % Calculate J_loss
                    Jloss = J_loss_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting);

                    % Calculate the loss angle and save these results to the
                    % structure variable
                    Jabs = sqrt(J_storage_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting).^2 + ...
                        J_loss_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting).^2);
                    theta = atan(((Jloss./(Jabs.^2))./(Jstorage./(Jabs.^2)))).*(180/pi);

                    dataStruct(indShift+i_loop).Jloss_conv = Jloss;
                    dataStruct(indShift+i_loop).Jstorage_conv = Jstorage;
                    dataStruct(indShift+i_loop).theta_conv = theta;

                    % Lastly, the Linear J parameters
                    % Calculate J_storage
                    Jstorage_linear = J_storage_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting);

                    % Calculate J_loss
                    Jloss_linear = J_loss_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting);

                    % Finally, calculate the loss angle and save these results to the
                    % structure variable
                    Jabs_linear = sqrt(J_storage_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting).^2 + ...
                        J_loss_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting).^2);
                    theta_linear = atan(((Jloss_linear./(Jabs_linear.^2))./(Jstorage_linear./(Jabs_linear.^2)))).*(180/pi);


                    dataStruct(indShift+i_loop).Jloss_linear = Jloss_linear;
                    dataStruct(indShift+i_loop).Jstorage_linear = Jstorage_linear;
                    dataStruct(indShift+i_loop).theta_linear = theta_linear;

                case 'open'

                    if timeOperation
                        iterativeFitStartTime(iVoigt,i_loop) = toc;
                    end

                    % Make the Generalized Voigt for this loop
                    [U_func,F_conv,F_conv_wrapper,lb,ub,subref,selector] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);

                    % Make Schapery Functions
    %                 [~,E_star,A_conv,F_conv_schapery_area,F_conv_schapery_area_wrapper] = makeSchaperyModel(elasticSetting,fluidSetting,'n',timeVec,dt,st,minTimescale,'voigt',n_terms,advancedLog,forwardFitTimescale,r_tip,E_tip,nu_tip,nu_sample);
                    [~,E_star,A_conv,h_conv_schapery_area,h_conv_schapery_area_wrapper] = makeSchaperyModel(elasticSetting_maxwell,fluidSetting_maxwell,'n',maxwellTimeVec,dt,st,minTimescale,'maxwell',n_loss_fit,advancedLog,forwardFitTimescale,r_tip,E_tip,nu_tip,nu_sample);

                    % Function Shape for Multiple Voigt with Linear Load Assumption
                    [U_linear,F_linear_fit] = makeLinearViscoFunction(elasticSetting,fluidSetting,n_terms,timeVec,dt,st);

                    % Create the Gen. Maxwell Force Convolution Model
                    [G_func_maxwell,Q_func_maxwell,h_conv_maxwell,h_conv_maxwell_wrapper,padSize] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_loss_fit,maxwellTimeVec,dt,st);

                    % Create the Harmonic Gen. Maxwell Models
                    [E_storage,E_loss,lossAngle,harmonic_wrapper,lb_maxwell,ub_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_loss_fit,maxwellTimeVec,minTimescale,padSize,forwardFitTimescale);

                    [tauInds,modulusInds] = getParamIndices(ub_loop,elasticSetting,fluidSetting);
                    [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(ub_loop_maxwell,elasticSetting_maxwell,fluidSetting_maxwell);

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
                    
                    beta_dist = zeros(length(ub),n_samples);
                    beta0_dist = zeros(size(beta_dist));
                    resnorm_dist = zeros(1,n_samples);
                    beta_dist_maxwell = zeros(length(ub_maxwell),n_samples);
                    beta0_dist_maxwell = zeros(size(beta_dist_maxwell));
                    resnorm_dist_maxwell = zeros(1,n_samples);
                    beta_dist_schapery_area = zeros(length(ub_maxwell),n_samples);
                    beta0_dist_schapery_area = zeros(size(beta_dist_schapery_area));
                    resnorm_dist_schapery_area = zeros(1,n_samples);
                    beta_dist_linear = zeros(length(ub),n_samples);
                    beta0_dist_linear = zeros(size(beta_dist_linear));
                    resnorm_dist_linear = zeros(1,n_samples);

                    y_fit = dataStruct(indShift+i_loop).h_norm_log;
                    y_fit_maxwell = dataStruct(indShift+i_loop).F_norm_log;
                    y_fit_linear = dataStruct(indShift+i_loop).chi_linear_log;
                    y_fit_area = log_scale(F_schapery_area,timeVec,dt,st);

                    if timeOperation
                        parallelFitStartTime(iVoigt,i_loop,n_terms) = toc;
                        convTime = Par(n_samples);
                        maxwellTime = Par(n_samples);
                        schaperyTime = Par(n_samples);
                        if i == 1
                           convTimeTotal = [];
                           maxwellTimeTotal = [];
                           schaperyTimeTotal = [];
                        end
                        ticBytes(gcp);
                        if ~exist('userMem','var')
                            %[userMem,systemMem] = memory;
                            userMem = 0;
                            systemMem = 0;
                        end
                        ramused = NaN(1,n_samples);
                        ramav = NaN(1,n_samples);
                        ramtime = NaN(1,n_samples);
                    end

                    manualVoigtInputParam = 0;
                    manualMaxwellInputParam = 0;

                    if weightedFit
                        
                        % New weighting based on dataset
%                         weightVectorVoigt = (1./h_norm_log_loop);
%                         weightVectorMaxwell = (1./F_norm_log_loop);
%                         weightVectorSchapery = (1./F_log_loop);
%                         weightVectorDataVoigt = (1./h_norm_log_loop).*(1./(1-exp(t_log_loop.*((nu_sample-1)./(dt)))));
%                         weightVectorDataMaxwell = (1./F_norm_log_loop).*((1-exp(t_log_loop.*((nu_sample-1)./(dt)))));
%                         weightVectorDataSchapery = (1./F_log_loop).*(1+nu_sample).*((1-exp(t_log_loop.*((nu_sample-1)./(dt)))));
        
                        % time weighting
                        timeweight = (t_plot_loop./(max(t_plot_loop)*0.9));
                        timeweight(timeweight > 1) = 1;
                        weightVectorVoigt = timeweight;
                        weightVectorMaxwell = timeweight;
                        weightVectorSchapery = timeweight;
                        weightVectorDataVoigt = timeweight;
                        weightVectorDataMaxwell = timeweight;
                        weightVectorDataSchapery = timeweight;

                    else
                        weightVectorVoigt = 1;
                        weightVectorMaxwell = 1;
                        weightVectorSchapery = 1;
                        weightVectorDataVoigt = 1;
                        weightVectorDataMaxwell = 1;
                        weightVectorDataSchapery = 1;
                    end

                    % Parallel Loop for Random Search
                    parfor j = 1:n_samples
                        warning('off');
                        lastwarn(''); % Clear last warning

                        try
                            % Initialize
                            beta0_linear = [];
                            beta0_conv = [];
                            beta0_maxwell = [];
                            beta0_schapery = [];

                            beta0_linear = NaN(size(ub));
                            beta0_maxwell = NaN(size(ub_maxwell));
                            beta0_conv = NaN(size(ub));
                            beta0_schapery = NaN(size(ub));

                            ub_randLimit = 0;
                            lb_randLimit = -6;
                            ub_randLimit_maxwell = 6;
                            lb_randLimit_maxwell = 0;
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

                            opts = optimoptions('fmincon','Algorithm','active-set',...
                                'MaxFunctionEvaluations', n_maxIterations,...
                                'MaxIterations', n_maxIterations,...
                                'Display','none',...
                                'FiniteDifferenceType','central',...
                                'FunctionTolerance', 0,...
                                'OptimalityTolerance', 0,...
                                'StepTolerance', 0);

                            % For the Voigt models
                            beta0_loop = ((ub-lb).*rand(size(ub)) + lb);
                            if elasticSetting == 'y'
                                beta0_loop(1) = 10^((ub_randLimit-lb_randLimit)*rand() + lb_randLimit);
                                beta0_loop(2:2:end) = 10^((ub_randLimit-lb_randLimit)*rand() + lb_randLimit);
                            else
                                beta0_loop(1:2:end) = 10^((ub_randLimit-lb_randLimit)*rand(size(beta0_loop(1:2:end))) + lb_randLimit);
                            end
                            if fluidSetting == 'y'
                                beta0_loop(end) = 10^(ub_fluidityRandLimit_init*rand());
                            end

                            beta0_linear = beta0_loop;
                            beta0_conv = beta0_loop;
    %                         beta0_schapery = beta0_loop;

                            % For testing against ELopezSim
                            if manualVoigtInputParam
                                if strcmp(elasticSetting,'y')
                                    if strcmp(fluidSetting,'y')
                                        beta0_conv(1) = manualParams(1);
                                        beta0_conv(2:2:end-1) = manualParams(2:2:length(beta0_loop)-1);
                                        beta0_conv(3:2:end-1) = manualParams(3:2:length(beta0_loop)-1);
                                        beta0_conv(end) = manualParams(end);
                                    else
                                        beta0_conv(1) = manualParams(1);
                                        beta0_conv(2:2:end) = manualParams(2:2:length(beta0_loop));
                                        beta0_conv(3:2:end) = manualParams(3:2:length(beta0_loop));
                                    end
                                else
                                    if strcmp(fluidSetting,'y')
                                        beta0_conv(1:2:end-1) = manualParams(1:2:length(beta0_loop)-1);
                                        beta0_conv(2:2:end-1) = manualParams(2:2:length(beta0_loop)-1);
                                        beta0_conv(end) = manualParams(end);
                                    else
                                        beta0_conv(1:2:end) = manualParams(1:2:length(beta0_loop));
                                        beta0_conv(2:2:end) = manualParams(2:2:length(beta0_loop));
                                    end
                                end
                            end

                            % Re-do for the maxwell/schapery
                            beta0_loop = ((ub_maxwell-lb_maxwell).*rand(size(ub_maxwell)) + lb_maxwell);
                            if elasticSetting_maxwell == 'y'
                                beta0_loop(1) = 10^((ub_randLimit_maxwell-lb_randLimit_maxwell)*rand() + lb_randLimit_maxwell);
                                beta0_loop(2:2:end) = 10^((ub_randLimit_maxwell-lb_randLimit_maxwell)*rand() + lb_randLimit_maxwell);
                            else
                                beta0_loop(1:2:end) = 10^((ub_randLimit_maxwell-lb_randLimit_maxwell)*rand(size(beta0_loop(1:2:end))) + lb_randLimit_maxwell);
                            end
                            if fluidSetting_maxwell == 'y'
                                beta0_loop(end) = 10^(ub_fluidityRandLimit_init*rand());
                            end

                            beta0_maxwell = beta0_loop;
                            beta0_schapery = beta0_loop;

                            if timeOperation
                                Par.tic;
                            end

                            if ~fitWithLSQ

                                % Fit the GKV model
    %                             Fsumsquares = @(c) sum(((F_conv_wrapper(c,x_fit) - abs(y_fit))./abs(y_fit)).^2);
    %                             Fsumsquares = @(c) sum((F_conv_wrapper(c,x_fit) - y_fit).^2);
                                Fsumsquares = @(c) sum((weightVectorVoigt.*F_conv_wrapper(c,x_fit) - weightVectorDataVoigt.*y_fit).^2);

                                opts.StepTolerance = scaling;
    %                             opts.Display = 'none';
    %                             opts.OutputFcn = @outfun;

                                [betatemp,ressquared,eflag,outputu] = ...
                                    fmincon(Fsumsquares,beta0_conv,[],[],[],[],lb,ub,[],opts);

                                opts.StepTolerance = 0;
    %                             opts.Display = 'none';
    %                             opts.OutputFcn = [];

                                residual = (F_conv_wrapper(betatemp,x_fit) - y_fit);

                            else
                                lsqoptions.StepTolerance = scaling;
    %                             lsqoptions.Display = 'none';
    %                             lsqoptions.MaxFunctionEvaluations = 1e3;
    %                             lsqoptions.OutputFcn = @outfun;

                                if ~fitSecondDatastream 
                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                        lsqcurvefit(@(c,inputs) weightVectorVoigt.*F_conv_wrapper(c,inputs),beta0_conv,x_fit,weightVectorDataVoigt.*y_fit,lb,ub,lsqoptions);
                                else
                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                        lsqcurvefit(voigtFitWrapper,beta0_conv,x_fit,y_fit_alt,lb,ub,lsqoptions);
                                end
                                lsqoptions.StepTolerance = 0;
    %                             lsqoptions.Display = 'none';
    %                             lsqoptions.MaxFunctionEvaluations = n_maxIterations;
    %                             lsqoptions.OutputFcn = [];
                            end

                            voigtData = [log_scale(x_fit(:,1),x_fit(:,1),dt,st)',F_conv_wrapper(betatemp,x_fit)'];

                            if timeOperation
                                convTime(j) = Par.toc;
                            end

                            if ~ignoreWarns
                                [warnMsg, warnId] = lastwarn;
                                if ~isempty(warnMsg)
                                    if contains(warnMsg,'inaccurate')
                                        beta_dist(:,j) = NaN;
                                        beta0_dist(:,j) = NaN;
                                        resnorm_dist(j) = NaN;

                                        beta_dist_schapery_area(:,j) = NaN;
                                        beta0_dist_schapery_area(:,j) = NaN;
                                        resnorm_dist_schapery_area(j) = NaN;

                                        beta_dist_linear(:,j) = NaN;
                                        beta0_dist_linear(:,j) = NaN;
                                        resnorm_dist_linear(j) = NaN;

                                        fprintf('Entered a NaN Row for ALL cases.\n');
                                        continue;
                                    end
                                end
                            end

                            if ~fitSecondDatastream 
    %                             resnorm_dist(j) = sum((residual./y_fit).^2);
    %                             resnorm_dist(j) = sum((residual).^2);
                                resnorm_dist(j) = sum(((F_conv_wrapper(betatemp,x_fit)-y_fit).^2)./movvar(y_fit,3))./(length(y_fit)-length(betatemp)); % Standard Error of Regression (S)
    %                             resnorm_dist(j) = sum(abs(((F_conv_wrapper(betatemp,x_fit)-y_fit))./y_fit));
                            else
                                resnorm_dist(j) = sum((residual./y_fit_alt).^2,'all');
                            end

                            if timeOperation
                                Par.tic;
                            end

                            if ~fitWithLSQ

                                % Fit the Maxwell model
    %                             Fsumsquares_maxwell = @(c) sum(abs(h_conv_maxwell_wrapper(c,maxwellInputs) - y_fit_maxwell)./y_fit_maxwell);
    %                             Fsumsquares_maxwell = @(c) sum((h_conv_maxwell_wrapper(c,maxwellInputs) - y_fit_maxwell).^2);
                                Fsumsquares_maxwell = @(c) sum((weightVectorMaxwell.*h_conv_maxwell_wrapper(c,maxwellInputs) - weightVectorDataMaxwell.*y_fit_maxwell).^2);

                                opts.StepTolerance = scaling;
    %                             opts.Display = 'none';
    %                             opts.OutputFcn = @outfun;

                                [betatemp,ressquared,eflag,outputu] = ...
                                    fmincon(Fsumsquares_maxwell,beta0_maxwell,[],[],[],[],lb_maxwell,ub_maxwell,[],opts);

                                opts.StepTolerance = 0;
    %                             opts.Display = 'none';
    %                             opts.OutputFcn = [];

                                residual = (h_conv_maxwell_wrapper(betatemp,maxwellInputs) - y_fit_maxwell);

                            else

                                lsqoptions.StepTolerance = scaling;
    %                             lsqoptions.Display = 'none';
    %                             lsqoptions.MaxFunctionEvaluations = 1e3;
    %                             lsqoptions.OutputFcn = @outfun;

                                if ~fitSecondDatastream 
    %                                     [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    %                                         lsqcurvefit(@(c,inputs) 0.1*weightedFit+weightVectorMaxwell.*h_conv_maxwell_wrapper(c,inputs),beta0_maxwell,maxwellInputs,weightVectorDataMaxwell.*y_fit_maxwell,lb_maxwell,ub_maxwell,lsqoptions);
                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                        lsqcurvefit(@(c,inputs) 0.1*weightedFit+weightVectorMaxwell.*h_conv_maxwell_wrapper(c,inputs),beta0_maxwell,maxwellInputs,weightVectorDataMaxwell.*y_fit_maxwell,lb_maxwell,ub_maxwell,lsqoptions);
                                else
                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                        lsqcurvefit(maxwellFitWrapper,beta0_maxwell,maxwellInputs,y_fit_maxwell_alt,lb_maxwell,ub_maxwell,lsqoptions);
                                end

                                if exitflag == 1
                                    disp('Maxwell Model Converged!');   
                                end

                                lsqoptions.StepTolerance = 0;
    %                             lsqoptions.Display = 'none';
    %                             lsqoptions.MaxFunctionEvaluations = n_maxIterations;
    %                             lsqoptions.OutputFcn = [];

                            end

                            maxwellData = [log_scale(maxwellInputs(:,1),maxwellInputs(:,1),dt,st)',h_conv_maxwell_wrapper(betatemp,maxwellInputs)'];

                            if timeOperation
                                maxwellTime(j) = Par.toc;
                            end

                            if ~ignoreWarns
                                [warnMsg, warnId] = lastwarn;
                                if ~isempty(warnMsg)
                                    if contains(warnMsg,'inaccurate')
                                        beta_dist(:,j) = NaN;
                                        beta0_dist(:,j) = NaN;
                                        resnorm_dist(j) = NaN;

                                        fprintf('Entered a NaN Row.\n');
                                        continue;
                                    end
                                end
                            end

                            beta_dist_maxwell(:,j) = betatemp';
                            beta0_dist_maxwell(:,j) = beta0_maxwell;
                            if ~fitSecondDatastream 
    %                             resnorm_dist_maxwell(j) = sum((residual./y_fit_maxwell).^2);
    %                             resnorm_dist_maxwell(j) = sum((residual).^2);
                                resnorm_dist_maxwell(j) = sum(((h_conv_maxwell_wrapper(betatemp,maxwellInputs)-y_fit_maxwell).^2)./movvar(y_fit_maxwell,3))./(length(y_fit_maxwell)-length(betatemp)); % Standard Error of Regression (S)
    %                             resnorm_dist_maxwell(j) = sum(abs(((h_conv_maxwell_wrapper(betatemp,maxwellInputs)-y_fit_maxwell))./y_fit_maxwell));
                            else
                                resnorm_dist_maxwell(j) = sum((residual./y_fit_maxwell_alt).^2,'all');
                            end

                            if timeOperation
                                Par.tic;
                            end

                            if ~fitWithLSQ

                                % Fit the Maxwell model
    %                             Fsumsquaresschapery = @(c) sum(abs(h_conv_schapery_area_wrapper(c,x_fit_area) - y_fit_area)./y_fit_area);
    %                             Fsumsquaresschapery = @(c) sum((h_conv_schapery_area_wrapper(c,x_fit_area) - y_fit_area).^2);
                                Fsumsquaresschapery = @(c) sum((weightVectorSchapery.*h_conv_schapery_area_wrapper(c,x_fit_area) - weightVectorDataSchapery.*y_fit_area).^2);

                                opts.StepTolerance = scaling;
    %                             opts.Display = 'none';
    %                             opts.OutputFcn = @outfun;

                                [betatemp,ressquared,eflag,outputu] = ...
                                    fmincon(Fsumsquaresschapery,beta0_schapery,[],[],[],[],lb_maxwell,ub_maxwell,[],opts);

                                opts.StepTolerance = 0;
    %                             opts.Display = 'none';
    %                             opts.OutputFcn = [];

                                residual = (h_conv_schapery_area_wrapper(betatemp,x_fit_area) - y_fit_area);

                            else

                                lsqoptions.StepTolerance = scaling;
    %                                 lsqoptions.MaxFunctionEvaluations = 1e3;
    %                                 lsqoptions.Display = 'none';
    %                                 lsqoptions.OutputFcn = @outfun;

                                if ~fitSecondDatastream
                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                        lsqcurvefit(@(c,inputs) weightVectorSchapery.*h_conv_schapery_area_wrapper(c,inputs),beta0_schapery,x_fit_area,weightVectorDataSchapery.*y_fit_area,lb_maxwell,ub_maxwell,lsqoptions);
                                else
                                    [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                                        lsqcurvefit(schaperyFitWrapper,beta0_schapery,x_fit_area,y_fit_area_alt,lb_maxwell,ub_maxwell,lsqoptions);
                                end

                                if exitflag == 1
                                    disp('Schapery-Maxwell Model Converged!');
                                end

                                lsqoptions.StepTolerance = 0;
    %                             lsqoptions.MaxFunctionEvaluations = n_maxIterations;
    %                             lsqoptions.Display = 'none';
    %                             lsqoptions.OutputFcn = @outfun;

                            end

                            schaperyData = [log_scale(x_fit_area(:,1),x_fit_area(:,1),dt,st)',h_conv_schapery_area_wrapper(betatemp,x_fit_area)'];

                            if timeOperation
                                schaperyTime(j) = Par.toc;
                            end

                            if ~ignoreWarns
                                [warnMsg, warnId] = lastwarn;
                                if ~isempty(warnMsg)
                                    if contains(warnMsg,'inaccurate')
                                        beta_dist_schapery_area(:,j) = NaN;
                                        beta0_dist_schapery_area(:,j) = NaN;
                                        resnorm_dist_schapery_area(j) = NaN;

                                        beta_dist_linear(:,j) = NaN;
                                        beta0_dist_linear(:,j) = NaN;
                                        resnorm_dist_linear(j) = NaN;

                                        fprintf('Entered a NaN Row for the Schapery and Linear cases.\n');
                                        continue;
                                    end
                                end
                            end

                            beta_dist_schapery_area(:,j) = betatemp';
                            beta0_dist_schapery_area(:,j) = beta0_schapery;
                            if ~fitSecondDatastream
    %                             resnorm_dist_schapery_area(j) = sum((residual./y_fit_area).^2);
    %                             resnorm_dist_schapery_area(j) = sum((residual).^2);
                                resnorm_dist_schapery_area(j) = sum(((h_conv_schapery_area_wrapper(betatemp,x_fit_area)-y_fit_area).^2)./movvar(y_fit_area,3))./(length(y_fit_area)-length(betatemp)); % Standard Error of Regression (S)
    %                             resnorm_dist_schapery_area(j) = sum(abs(((h_conv_schapery_area_wrapper(betatemp,x_fit_area)-y_fit_area))./y_fit_area));
                            else
                                resnorm_dist_schapery_area(j) = sum((residual./y_fit_area_alt).^2,'all');
                            end

                            % Fit the linear model
    %                         Fsumsquares_linear = @(c) sum(((F_linear_fit(c,timeVec) - abs(y_fit_linear))./abs(y_fit_linear)));
                            Fsumsquares_linear = @(c) sum(abs(F_linear_fit(c,timeVec) - y_fit_linear));
                            opts = optimoptions('fmincon','Algorithm','interior-point',...
                                'Display','none','FiniteDifferenceType','central',...
                                'FunctionTolerance', scaling,...
                                'OptimalityTolerance',scaling,...
                                'StepTolerance', scaling);

                            [betatemp,ressquared,eflag,outputu] = ...
                                fmincon(Fsumsquares_linear,beta0_linear,[],[],[],[],lb,ub,[],opts);

                            if ~ignoreWarns
                                [warnMsg, warnId] = lastwarn;
                                if ~isempty(warnMsg)
                                    if contains(warnMsg,'inaccurate')
                                        beta_dist_linear(:,j) = NaN;
                                        beta0_dist_linear(:,j) = NaN;
                                        resnorm_dist_linear(j) = NaN;

                                        fprintf('Entered a NaN Row for the Linear case.\n');
                                        continue;
                                    end
                                end
                            end

                            beta_dist_linear(:,j) = betatemp';
                            beta0_dist_linear(:,j) = beta0_linear;
                            resnorm_dist_linear(j) = ressquared;

                            if timeOperation
                                % Track memory
                                %ramused(j) = userMem.MemUsedMATLAB;
                                ramused(j) = 0;
                                %ramav(j) = systemMem.PhysicalMemory.Available;
                                ramav(j) = 0;
                                ramtime(j) = now;
                            end

                        catch ERROR
                            beta_dist(:,j) = NaN;
                            beta0_dist(:,j) = NaN;
                            resnorm_dist(j) = NaN;

                            beta_dist_maxwell(:,j) = NaN;
                            beta0_dist_maxwell(:,j) = NaN;
                            resnorm_dist_maxwell(j) = NaN;

                            beta_dist_schapery_area(:,j) = NaN;
                            beta0_dist_schapery_area(:,j) = NaN;
                            resnorm_dist_schapery_area(j) = NaN;

                            beta_dist_linear(:,j) = NaN;
                            beta0_dist_linear(:,j) = NaN;
                            resnorm_dist_linear(j) = NaN;

                            if timeOperation
                                convTime(j) = 1;
                                maxwellTime(j) = 1;
                                schaperyTime(j) = 1;
                            end

                            voigtData = [0,0];
                            maxwellData = [0,0];
                            schaperyData = [0,0];

                            fprintf('\nTerm Number %d Search: Entered a NaN Row for loop #%d because there was an ERROR:\n',n_terms,j);
                            disp(ERROR.message)
                        end

                    end

                    if timeOperation
                        parallelFitEndTime(iVoigt,i_loop,n_terms) = toc;
                        parallelFitParallelBytes(iVoigt,i_loop,n_terms) = getfield(tocBytes(gcp),{1});

                        stop(convTime);
                        stop(maxwellTime);
                        stop(schaperyTime);

                        convTimeTotal = convTime;
                        maxwellTimeTotal = maxwellTime;
                        schaperyTimeTotal = schaperyTime;

                        fprintf('\nThe parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                            (parallelFitEndTime(iVoigt,i_loop,n_terms)-parallelFitStartTime(iVoigt,i_loop,n_terms))/60, parallelFitParallelBytes(iVoigt,i_loop,n_terms)*1e-6);
                    end

                    n_ranks = 10; % Number of places to save the parameter sets of

    %                 if n_ranks < length(beta_dist(1,~isnan(beta_dist(1,:))))
    %                     n_ranks = length(beta_dist(1,~isnan(beta_dist(1,:))));
    %                 end

                    resnorm_dist_linear_temp = resnorm_dist_linear;
                    resnorm_dist_linear_temp(resnorm_dist_linear_temp == 0) = NaN;
                    [~,idx_linear] = min(resnorm_dist_linear_temp,[],'omitnan');
                    dataStruct(indShift+i_loop).linearParams = beta_dist_linear(:,idx_linear);

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

                    resnorm_dist_schapery_area_temp = resnorm_dist_schapery_area;
                    resnorm_dist_schapery_area_temp(resnorm_dist_schapery_area_temp == 0) = NaN;
                    [~,idx_schapery] = min(resnorm_dist_schapery_area_temp,[],'omitnan');
                    dataStruct(indShift+i_loop).schaperyParams = beta_dist_schapery_area(:,idx_schapery);

                    [~,I_schapery] = sort(resnorm_dist_schapery_area_temp,'MissingPlacement','last');
                    dataStruct(indShift+i_loop).schaperyParams_alternatives = beta_dist_schapery_area(:,I_schapery(1:n_ranks)); 

                    if timeOperation
                        iterativeFitEndTime(iVoigt,i_loop) = toc;
                        parallelFitBreakdown{iVoigt,i_loop} = {convTimeTotal,maxwellTimeTotal,schaperyTimeTotal};
                        parallelFitMemoryMonitor{iVoigt,i_loop} = vertcat(ramused,ramav,ramtime);
                    end

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
                    if elasticSetting == 'y'
                        fprintf(fid,'E_e = %3.6g [Pa]\r\n',beta_dist_maxwell(1,idx_maxwell));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,idx_maxwell),iii,beta_dist_maxwell(iii*2+1,idx_maxwell));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,idx_maxwell),iii,beta_dist_maxwell(iii*2,idx_maxwell));
                        end
                        if fluidSetting == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

    %                 fprintf(fid,'\r\nSchapery Generalized Kelvin-Voigt\n');
    %                 if elasticSetting == 'y'
    %                     fprintf(fid,'J_g = %3.6g [1/Pa]\r\n',beta_dist_schapery_area(1,idx_schapery));
    % 
    %                     for iii = 1:n_terms
    %                         fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,idx_schapery),iii,beta_dist_schapery_area(iii*2+1,idx_schapery));
    %                     end
    %                     if fluidSetting == 'y'
    %                         fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
    %                     else
    %                         fprintf(fid,'No fluidity response used.\r\n');
    %                     end
    % 
    %                 else
    %                     fprintf(fid,'No elastic response used.\r\n');
    % 
    %                     for iii = 1:n_terms
    %                         fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,idx_schapery),iii,beta_dist_schapery_area(iii*2,idx_schapery));
    %                     end
    %                     if fluidSetting == 'y'
    %                         fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
    %                     else
    %                         fprintf(fid,'No fluidity response used.\r\n');
    %                     end
    % 
    %                 end

                    fprintf(fid,'\r\n\r\nSchapery Generalized Maxwell\r\n');
                    if elasticSetting_maxwell == 'y'
                        fprintf(fid,'E_e = %3.6g [Pa]\r\n',beta_dist_schapery_area(1,idx_schapery));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,idx_schapery),iii,beta_dist_schapery_area(iii*2+1,idx_schapery));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,idx_schapery),iii,beta_dist_schapery_area(iii*2,idx_schapery));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
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
                            fprintf(fid,'E_g = %3.6g [Pa]\r\n',beta_dist_maxwell(1,I_maxwell(alti)));

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2+1,I_maxwell(alti)));
                            end
                            if fluidSetting == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        else
                            fprintf(fid,'No elastic response used.\r\n');

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2,I_maxwell(alti)));
                            end
                            if fluidSetting == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        end

    %                     % Print Schapery
    %                     fprintf(fid,'\r\nSchapery Generalized Kelvin-Voigt\n');
    %                     if elasticSetting == 'y'
    %                         fprintf(fid,'J_g = %3.6g [1/Pa]\r\n',beta_dist_schapery_area(1,I_schapery(alti)));
    % 
    %                         for iii = 1:n_terms
    %                             fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2+1,I_schapery(alti)));
    %                         end
    %                         if fluidSetting == 'y'
    %                             fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
    %                         else
    %                             fprintf(fid,'No fluidity response used.\r\n');
    %                         end
    % 
    %                     else
    %                         fprintf(fid,'No elastic response used.\r\n');
    % 
    %                         for iii = 1:n_terms
    %                             fprintf(fid,'J_%d = %3.6g [1/Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2,I_schapery(alti)));
    %                         end
    %                         if fluidSetting == 'y'
    %                             fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
    %                         else
    %                             fprintf(fid,'No fluidity response used.\r\n');
    %                         end
    % 
    %                     end                    

                        % Print Schapery Maxwell
                        fprintf(fid,'\r\nSchapery Generalized Maxwell\n');
                        if elasticSetting_maxwell == 'y'
                            fprintf(fid,'E_g = %3.6g [Pa]\r\n',beta_dist_schapery_area(1,I_schapery(alti)));

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2+1,I_schapery(alti)));
                            end
                            if fluidSetting_maxwell == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        else
                            fprintf(fid,'No elastic response used.\r\n');

                            for iii = 1:n_terms
                                fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2,I_schapery(alti)));
                            end
                            if fluidSetting_maxwell == 'y'
                                fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
                            else
                                fprintf(fid,'No fluidity response used.\r\n');
                            end

                        end

                    end

                    % Calculate Harmonic Quantities
                    if smoothData == 'n'
                        inputs = [dataStruct(indShift+i_loop).t_r, dataStruct(indShift+i_loop).F_r];
                    else
                        inputs = [dataStruct(indShift+i_loop).t_r_smooth, dataStruct(indShift+i_loop).F_r_smooth];
                    end

    %                 dt = (dataStruct(indShift+i_loop).dt);
    %                 st = timeVec(end);
    %                 t_log = log_scale(timeVec,timeVec,dt,st);

                    % Define the convolution fit dataset generation function
                    data_fit = cumtrapz(F_conv(dataStruct(indShift+i_loop).convParams,inputs)).*dt;
                    data_fit_log = log_scale(data_fit, timeVec, dt, st);

                    % Define the linear convolution fit dataset generation function
                    linear_fit = cumtrapz(F_conv(dataStruct(indShift+i_loop).linearParams,inputs)).*dt;
                    linear_fit_log = log_scale(linear_fit, timeVec, dt, st);

                    % Move on to plotting the loss angle
                    % First, generate a frequency array in log scale
                    de0 = 2.0.*pi.*(1./timeVec(end));
                    maxi = 2.0.*pi.*(1./mean(diff(timeVec)));
                    omega = log_tw(de0,maxi);

                    % Second, use the convolution J parameters
                    % Calculate J_storage
                    Jstorage = J_storage_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting);

                    % Calculate J_loss
                    Jloss = J_loss_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting);

                    % Calculate the loss angle and save these results to the
                    % structure variable
                    Jabs = sqrt(J_storage_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting).^2 + ...
                        J_loss_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting).^2);
                    theta = atan(((Jloss./(Jabs.^2))./(Jstorage./(Jabs.^2)))).*(180/pi);

                    dataStruct(indShift+i_loop).Jloss_conv = Jloss;
                    dataStruct(indShift+i_loop).Jstorage_conv = Jstorage;
                    dataStruct(indShift+i_loop).theta_conv = theta;

                    % Lastly, the Linear J parameters
                    % Calculate J_storage
                    Jstorage_linear = J_storage_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting);

                    % Calculate J_loss
                    Jloss_linear = J_loss_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting);

                    % Finally, calculate the loss angle and save these results to the
                    % structure variable
                    Jabs_linear = sqrt(J_storage_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting).^2 + ...
                        J_loss_advanced(omega,dataStruct(indShift+i_loop).linearParams,elasticSetting,fluidSetting).^2);
                    theta_linear = atan(((Jloss_linear./(Jabs_linear.^2))./(Jstorage_linear./(Jabs_linear.^2)))).*(180/pi);

                    dataStruct(indShift+i_loop).Jloss_linear = Jloss_linear;
                    dataStruct(indShift+i_loop).Jstorage_linear = Jstorage_linear;
                    dataStruct(indShift+i_loop).theta_linear = theta_linear;

                otherwise
                    error('The fit method provided is not implemented. Please check the value of fitMethod.')

            end % End switch fitMethod

            if timeOperation
                fitSwitchEndTime(iVoigt,i_loop) = toc;
            end

            % Check the performance
            if plotFigs
                figure(1)
                clf
                scatter(t_log,y_fit,'r');
                hold on
                scatter(t_log,y_fit_area,'b');
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                title('Convolution-Method Force Fit')
                xlabel('Time [s]')
                ylabel('$\int_{0}^{t} U(t-\xi) F(\xi) d\xi$ [$m^2$]', 'Interpreter', 'latex', 'FontSize', 16)
                plot(t_log,F_conv_wrapper(beta_dist(:,idx),x_fit),'r');
                plot(t_log,h_conv_schapery_area_wrapper(beta_dist_schapery_area(:,idx_schapery),x_fit),'b');
                grid on
                tempString = sprintf('Convolution Fit\n%d-Terms, v_{approach} [nm/s] = %1.4g',n_terms,v_unique(i_loop));
                tempString2 = sprintf('Area-Normal Schapery Fit\n%d-Terms, v_{approach} [nm/s] = %1.4g',n_terms,v_unique(i_loop));
                legend({'Averaged Data','Area-Normalized Average Data',tempString,tempString2},'Location','southoutside',...
                    'Orientation','horizontal')
                hold off

                saveas(gcf, [pathname sprintf('/AveragedConvFit-nVoigt_%d-%s',voigtArray(iVoigt),saveLabel)], 'jpg')
            end

            % Plot the first figure using the AFM data, and both new datasets
            if plotFigs
                figure(2)
                clf
                plot(t_log, dataStruct(indShift+i_loop).h_norm_log, 'k*', 'LineWidth', 3)
                hold on
                grid on
                plot(t_log, data_fit_log, 'c-', 'LineWidth', 5)
                plot(t_log, linear_fit_log, 'r--', 'LineWidth', 3)
                legend({'Experimental Data' 'Non-linear squares fit, Non Linear' 'Non-linear squares fit, Linear'}, 'Location', 'eastoutside', 'Orientation', 'vertical')
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                title('Non-linear Squares Fit', 'FontSize', 16)
                xlabel('Time [s]', 'FontSize', 20)
                ylabel('$\int_{0}^{t} U(t-\xi) F(\xi) d\xi$ [$m^2$]', 'Interpreter', 'latex', 'FontSize', 16)
                hold off

                saveas(gcf, [pathname sprintf('/AveragedDataFit-nVoigt_%d-%s',voigtArray(iVoigt),saveLabel)], 'jpg')
            end

            if plotFigs
                figure(3)
                clf
                semilogx(omega, theta, 'c', 'LineWidth', 5)
                hold on
                grid on
                semilogx(omega, theta_linear, 'r--', 'LineWidth', 5)
                legend({'Fit Non Linear' 'Fit Linear'})
                xlabel('$\omega, \,rad/s$', 'FontSize', 16, 'Interpreter', 'Latex')
                ylabel('$\delta(\omega), \,deg$', 'FontSize', 16, 'Interpreter', 'Latex')
                xlim([10 10^4])
                hold off

                saveas(gcf, [pathname sprintf('/AveragedLossAngle-nVoigt_%d-%s',voigtArray(iVoigt),saveLabel)], 'jpg')
            end

            if plotFigs
                figure(4)
                clf
                semilogx(omega, dataStruct(indShift+i_loop).theta_conv, 'k', 'LineWidth', 3)
                hold on
                grid on
                semilogx(omega, lossAngle(omega,beta_dist_maxwell(:,idx_maxwell)), 'r-', 'LineWidth', 3)
                legendString = 'Voigt Model';
    %                 for iii = 1:size(dataStruct(indShift+i_loop).maxwellParams_alternatives,2)
    %                     semilogx(omega, lossAngle(omega,dataStruct(indShift+i_loop).maxwellParams_alternatives(:,iii)), 'r-', 'LineWidth', 3)
    %                 end
                legend({legendString 'Maxwell Model'})
                xlabel('$\omega, \,rad/s$', 'FontSize', 16, 'Interpreter', 'Latex')
                ylabel('$\delta(\omega), \,deg$', 'FontSize', 16, 'Interpreter', 'Latex')
                xlim([10 10^4])
                title(sprintf('Loss Angle Comparison for %d Elements, %s',voigtArray(iVoigt),saveLabel))
                hold off

                saveas(gcf, [pathname sprintf('/LossAngleFitting-nElements_%d-%s',voigtArray(iVoigt),saveLabel)], 'jpg')
                savefig(gcf, [pathname sprintf('/LossAngleFitting-nElements_%d-%s',voigtArray(iVoigt),saveLabel)],'compact')
            end

            % Create the last plot
            if plotFigs
                figure(5)
                clf
                scatter(t_log,dataStruct(indShift+i_loop).h_norm_log,'r');
                hold on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                title('Convolution-Method Force Fit')
                xlabel('Time [s]')
    %             ylabel('$\int_{0}^{t} U(t-\xi) F(\xi) d\xi$ [$m^2$]', 'Interpreter', 'latex', 'FontSize', 16)
                ylabel('Integration Value [Various]', 'Interpreter', 'latex', 'FontSize', 16)
                plot(t_log,F_conv_wrapper(dataStruct(indShift+i_loop).convParams,x_fit),'r');
                plot(t_log,h_conv_maxwell_wrapper(dataStruct(indShift+i_loop).maxwellParams,maxwellInputs),'k');
    %                 for iii = 1:size(dataStruct(indShift+i_loop).maxwellParams_alternatives,2)
    %                     plot(t_log,F_conv_maxwell_wrapper(dataStruct(indShift+i_loop).maxwellParams_alternatives(:,iii),maxwellInputs),'k');
    %                 end
                grid on
                tempString = sprintf('Generalized Voigt Fit\n%d-Terms, v_{approach} [nm/s] = %1.4g',n_terms,v_unique(i_loop));
                tempString2 = sprintf('Generalized Maxwell Fit\n%d-Terms, v_{approach} [nm/s] = %1.4g',n_loss_fit,v_unique(i_loop));
                legend({'Averaged Data',tempString,tempString2},'Location','southoutside',...
                    'Orientation','horizontal')
                hold off

                saveas(gcf, [pathname sprintf('/AveragedConvFit-nVoigt_%d-nMaxwell_%d-%s',voigtArray(iVoigt),n_loss_fit,saveLabel)], 'jpg')
            end

            % Plot the overall Loss Modulus and Storage Modulus for
            % comparison with Alex's NIH Data (Maxwell)
            if plotFigs
                figure(6)
                clf
                plot(omega,E_loss(omega,dataStruct(indShift+i_loop).maxwellParams),'r');
                hold on
                plot(omega,E_storage(omega,dataStruct(indShift+i_loop).maxwellParams),'b');
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                tempString = sprintf('Predicted Viscoelastic Moduli, Generalized Maxwell Fit %d-Terms, v_{approach} [nm/s] = %1.4g',n_loss_fit,v_unique(i_loop));
                title(tempString)
                xlabel('$\omega, \,rad/s$', 'FontSize', 16, 'Interpreter', 'Latex')
                ylabel('Modulus [Pa]')
                grid on
                legend({'Loss Modulus','Storage Modulus'},'Location','southoutside',...
                    'Orientation','horizontal')
                hold off

                saveas(gcf, [pathname sprintf('/StorageAndLossModuli-nVoigt_%d-nMaxwell_%d-%s',voigtArray(iVoigt),n_loss_fit,saveLabel)], 'jpg')
            end

            % Plot the overall Loss Modulus and Storage Modulus for
            % comparison with Alex's NIH Data (Voigt)
            if plotFigs
                Jabs = sqrt(J_storage_advanced(omega,simParams,simElasticSetting,simFluidSetting).^2 + ...
                    J_loss_advanced(omega,simParams,simElasticSetting,simFluidSetting).^2);
                figure(7)
                clf
                plot(omega,J_loss_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting)./(Jabs.^2),'r');
                hold on
                plot(omega,J_storage_advanced(omega,dataStruct(indShift+i_loop).convParams,elasticSetting,fluidSetting)./(Jabs.^2),'b');
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                tempString = sprintf('Predicted Viscoelastic Moduli, Generalized Voigt Fit %d-Terms, v_{approach} [nm/s] = %1.4g',voigtArray(iVoigt),v_unique(i_loop));
                title(tempString)
                xlabel('$\omega, \,rad/s$', 'FontSize', 16, 'Interpreter', 'Latex')
                ylabel('Modulus [Pa]')
                grid on
                legend({'Loss Modulus','Storage Modulus'},'Location','southoutside',...
                    'Orientation','horizontal')
                hold off

                saveas(gcf, [pathname sprintf('/StorageAndLossModuli-nVoigt_%d-%s',voigtArray(iVoigt),saveLabel)], 'jpg')
            end

            if strcmp(avAnalysis,'i')
                save([pathname sprintf('/%s-dataStruct-%d_terms-%s.mat',savePrepend,n_terms,saveLabel)],'dataStruct')
                fprintf('\nSaved the dataStruct to an m-file for %d terms on File No. %d\n',n_terms,i_loop);
            end

            if saveGDdata && strcmp(optimizationMethod,'gd')
                % Set plotting threshold
                plotThreshold = Inf; % Standard Error (S)
                bigFontSize = 28;
                mediumFontSize = 16;
                markersize = 100;
                nPlotsPerRow = 3;
                nRow = ceil(size(beta_dist,1)/nPlotsPerRow);

                for ij = 1:3

                    X_beta = [];
                    Cost_beta = [];
                    alphaBeta = [];
                    modelLabel = [];

                    switch ij
                        case 1
                            X_beta = (beta_dist(:,resnorm_dist<plotThreshold)');
                            Cost_beta = repmat(resnorm_dist(resnorm_dist<plotThreshold)',[1 size(beta_dist,1)]).*(1./alphaModel);
                            alphaBeta = (min(resnorm_dist(resnorm_dist<plotThreshold))./resnorm_dist(resnorm_dist<plotThreshold))';
                            modelLabel = 'GKV';

                            % Make labels
                            [tauInds,modulusInds] = getParamIndices(X_beta,elasticSetting,fluidSetting);
                            
                            symLabels = cell(size(X_beta,1),1);
                            if strcmp(fluidSetting,'y')
                                symLabels(tauInds) = horzcat(sprintfc('\\tau_{%d}',(1:(length(tauInds)-1))),sprintf('\\phi_{f}'));
                            else
                                symLabels(tauInds) = sprintfc('\\tau_{%d}',(1:(length(tauInds))));
                            end
                            if strcmp(elasticSetting,'y')
                                symLabels(modulusInds) = horzcat(sprintf('J_{g}'),sprintfc('J_{%d}',(1:(length(modulusInds)-1))));
                            else
                                symLabels(modulusInds) = {sprintfc('J_{%d}',(1:(length(modulusInds))))};
                            end

                            yLabelString = 'Standard Error, S [$m^{3/2}$]';

                        case 2
                            X_beta = (beta_dist_maxwell(:,resnorm_dist_maxwell<plotThreshold)');
                            Cost_beta = repmat(resnorm_dist_maxwell',[1 size(beta_dist_maxwell,1)]).*alphaModel;
                            alphaBeta = (min(resnorm_dist_maxwell(resnorm_dist_maxwell<plotThreshold))./resnorm_dist_maxwell(resnorm_dist_maxwell<plotThreshold))';
                            modelLabel = 'GM';

                            % Make labels
                            [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(X_beta,elasticSetting_maxwell,fluidSetting_maxwell);
                            
                            symLabels = cell(size(X_beta,1),1);
                            if strcmp(fluidSetting_maxwell,'y')
                                symLabels(tauInds) = horzcat(sprintfc('\\tau_{%d}',(1:(length(tauInds)-1))),sprintf('\\phi_{f}'));
                            else
                                symLabels(tauInds) = sprintfc('\\tau_{%d}',(1:(length(tauInds))));
                            end
                            if strcmp(elasticSetting_maxwell,'y')
                                symLabels(modulusInds) = horzcat(sprintf('E_{e}'),sprintfc('E_{%d}',(1:(length(modulusInds)-1))));
                            else
                                symLabels(modulusInds) = {sprintfc('E_{%d}',(1:(length(modulusInds))))};
                            end

                            yLabelString = 'Standard Error, S [$N$]';

                        case 3
                            X_beta = (beta_dist_schapery_area(:,resnorm_dist_schapery_area<plotThreshold)');
                            Cost_beta = repmat(resnorm_dist_schapery_area',[1 size(beta_dist_schapery_area,1)]).*alphaModel;
                            alphaBeta = (min(resnorm_dist_schapery_area(resnorm_dist_schapery_area<plotThreshold))./resnorm_dist_schapery_area(resnorm_dist_schapery_area<plotThreshold))';
                            modelLabel = 'SGM';

                            % Make labels
                            [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(X_beta,elasticSetting_maxwell,fluidSetting_maxwell);
                            
                            symLabels = cell(size(X_beta,1),1);
                            if strcmp(fluidSetting_maxwell,'y')
                                symLabels(tauInds) = horzcat(sprintfc('\\tau_{%d}',(1:(length(tauInds)-1))),sprintf('\\phi_{f}'));
                            else
                                symLabels(tauInds) = sprintfc('\\tau_{%d}',(1:(length(tauInds))));
                            end
                            if strcmp(elasticSetting_maxwell,'y')
                                symLabels(modulusInds) = horzcat(sprintf('E_{e}'),sprintfc('E_{%d}',(1:(length(modulusInds)-1))));
                            else
                                symLabels(modulusInds) = {sprintfc('E_{%d}',(1:(length(modulusInds))))};
                            end

                            yLabelString = 'Standard Error, S [$N$]';

                    end

                    % Begin visualizations of best parameter sets
                    if exist('paramSpacePlot','var')
                        clearvars paramSpacePlot
                    end
                    if exist('paramSpacePlot','var')
                        figure(paramSpacePlot);
                        clf;
                    else
                        paramSpacePlot = figure('position',[100 100 1200 500*nRow]);
                    end

                    for ii = 1:size(beta_dist,1)
                        if ii <= nPlotsPerRow
                            subtightplot (nRow, nPlotsPerRow, ii, [0.1 0.05], [0.2 0.1], [0.075 0.05]);
                        elseif ii > nPlotsPerRow && ii < nPlotsPerRow*(nRow)-(nPlotsPerRow-1)
                            subtightplot (nRow, nPlotsPerRow, ii, [0.1 0.05], [0.2 0.075], [0.075 0.05]);
                        else
                            subtightplot (nRow, nPlotsPerRow, ii, [0.1 0.05], [0.2 0.075], [0.075 0.05]);
                        end
                        hold on
                        for iii = 1:length(alphaBeta)
                            try
    %                                 scatter(Cost_beta(iii,ii),X_beta(iii,ii),markersize,'filled','MarkerFaceAlpha',alphaBeta(iii))
                                scatter(X_beta(iii,ii),Cost_beta(iii,ii),markersize,'filled','MarkerFaceAlpha',alphaBeta(iii))
                            catch ERROR
                            end
                        end
                        title(sprintf('$%s$',symLabels{ii}),'interpreter','latex', 'FontSize', bigFontSize)
                        if size(beta_dist,1)-ii < nPlotsPerRow
    %                             xlabel({'' '' 'Standard Error, S [$m^{3/2}$]'}, 'interpreter','latex', 'FontSize', mediumFontSize)
                            xlabel({'' '' 'Parameter Value [Various]'}, 'interpreter', 'latex', 'FontSize', mediumFontSize)
                        end
                        if ~mod(ii-1,nPlotsPerRow)
    %                                 ylabel('Parameter Value [Various]', 'interpreter', 'latex', 'FontSize', mediumFontSize)
                            ylabel(yLabelString, 'interpreter','latex', 'FontSize', mediumFontSize)
                        end
                        grid on
                        grid minor
                        set(gca,'XScale','log')
                        set(gca,'YScale','log')
                        hold off
                    end

                    saveas(paramSpacePlot, [pathname sprintf('/ParameterSpaceGD-nElements_%d-%s-%s',n_terms,modelLabel,saveLabel)], 'jpg')
                    savefig(paramSpacePlot, [pathname sprintf('/ParameterSpaceGD-nElements_%d-%s-%s',n_terms,modelLabel,saveLabel)],'compact')

                    close(paramSpacePlot);

                end

                if strcmp(avAnalysis,'a')
                    save([pathname sprintf('/%s-GDFittingInfo-%d_terms-LoadLevel%d-average.mat',savePrepend,n_terms,i_loop)],...
                        'beta0_dist','beta0_dist_maxwell','beta0_dist_schapery_area',...
                        'beta_dist','beta_dist_maxwell','beta_dist_schapery_area',...
                        'resnorm_dist','resnorm_dist_maxwell','resnorm_dist_schapery_area')
                    fprintf('\nSaved the GD Fitting Data to an m-file for %d terms and Load Level %d.\n',n_terms,i_loop);
                end

                if exist('paramSpacePlot','var')
                    clearvars paramSpacePlot
                end

            end

            fprintf(fid,'\r\n\r\n=================\r\n\r\n');

        end % End loopInd
    
    else
        
        % Make the save label
        switch avAnalysis
            case 'i'
                saveLabel = sprintf('AllFiles');
                endBound = length(Files);
                n_datasets = endBound;
            case 'a'
                saveLabel = sprintf('AllLoadLevels');
                endBound = size(dataStruct,2);
                n_datasets = endBound-(indShift);
        end

        x_lin = dataStruct(1+indShift).t_r_log; % X data for linear fit
        x_lin_mat = {dataStruct(1+indShift:endBound).t_r_log};
        x_lin_mat = cellfun(@transpose,x_lin_mat,'UniformOutput',false);
%         x_lin_mat = vertcat(x_lin_mat{:});
        y_lin = {dataStruct(1+indShift:endBound).F_r_log}; % Y data for linear fit
        y_lin_mat = {dataStruct(1+indShift:endBound).F_r_log};
        y_lin_mat = cellfun(@transpose,y_lin_mat,'UniformOutput',false);
        y_lin_mat = vertcat(y_lin_mat{:});

        init_params = polyfit(x_lin, y_lin{1}, 1); % Initial Linear Fit Parameters (guess)

        % Linear Fit Model Definition
        F_linFit = @(c,input) c .* input;
        F_linFit_str = '@(c,input) [';
        for j = 1:n_datasets
            F_linFit_str = horzcat(F_linFit_str,sprintf('c .* input{%d}',j));
            if j < n_datasets
                F_linFit_str = horzcat(F_linFit_str,';');
            else
                F_linFit_str = horzcat(F_linFit_str,']');
            end
        end
        F_linFit_mat = str2func(F_linFit_str);

        lb = 0; % Zero Newtons/s is the minimum force rate possible
        ub = 1e4; % 1e4 Newton/s is the upper bound

        options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
            'MaxFunctionEvaluations', n_maxIterations,...
            'MaxIterations', n_maxIterations,...
            'FiniteDifferenceType','central',...
            'FunctionTolerance', scaling,...
            'OptimalityTolerance', scaling,...
            'StepTolerance', scaling,...
            'Display', 'none');

        beta_dist = zeros(length(lb),n_samples);
        beta0_dist = zeros(size(beta_dist));
        resnorm_dist = zeros(1,n_samples);

        % Parallel Loop for Grid Search
        % Turned off parallel. To turn on, replace for with parfor.
        for i = 1:n_samples
            try
                beta0 = init_params(1);
                [betatemp,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                    lsqcurvefit(F_linFit_mat,beta0,x_lin_mat,y_lin_mat,lb,ub,options);

                beta_dist(:,i) = betatemp;
                beta0_dist(:,i) = beta0;
                resnorm_dist(i) = sum(abs(residual));
            catch
                beta_dist(:,i) = NaN;
                beta0_dist(:,i) = NaN;
                resnorm_dist(i) = NaN;
            end
            
        end

        [~,idx] = min(resnorm_dist);
        [dataStruct(1+indShift:endBound).Fdot] = deal(num2cell(repmat(beta_dist(1,idx),3,1)));

        clearvars x_lin x_lin_mat y_lin y_lin_mat beta_dist beta0_dist resnorm

        % Calculate the Linear Fluidity (chi)
        % According to Eq. 14, This is the relation between chi and tip
        % location when force is assumed to be linear
        dt = {(dataStruct(1+indShift:endBound).dt)};
        st = {};
        for i = 1:n_datasets
            st{i} = dataStruct(indShift+i).t_r(end);
        end

        n_loss_fit = n_terms;

        % Make the label for this analysis
        switch avAnalysis
            case 'i'
                fprintf(fid,'Simultaneous Fit for All Files:\r\n');
            case 'a'
                fprintf(fid,'Approach Velocities: %4.3g nm/s\r\n',v_unique);
        end

        if timeOperation
            timeStoreTemp = []; % Used to store the time for this loop
        end

        h_norm = {};
        F_norm = {};
        
        for i = 1:n_datasets
            if smoothData == 'n'
                h_norm{i} = ((8*sqrt(r_tip))/(3*(1-nu_sample))).*dataStruct(indShift+i).h_r.^1.5;
                F_norm{i} = ((3*(1-nu_sample))/(8*sqrt(r_tip))).*dataStruct(indShift+i).F_r;
            else
                h_norm{i} = ((8*sqrt(r_tip))/(3*(1-nu_sample))).*dataStruct(indShift+i).h_r_smooth.^1.5;
                F_norm{i} = ((3*(1-nu_sample))/(8*sqrt(r_tip))).*dataStruct(indShift+i).F_r_smooth;
            end
        end
        
        F_log = {};
        t_log = {};
        h_log = {};
        h_norm_log = {};

        % Get the input data organized
        if smoothData == 'n'
            for j = 1:n_datasets
                F_log{j} = log_scale(dataStruct(indShift+j).F_r,dataStruct(indShift+j).t_r,dt{j},st{j});
                t_log{j} = log_scale(dataStruct(indShift+j).t_r,dataStruct(indShift+j).t_r,dt{j},st{j});
                h_log{j} = log_scale(dataStruct(indShift+j).h_r,dataStruct(indShift+j).t_r,dt{j},st{j});
            
                F_schapery_area{j} = dataStruct(indShift+j).F_r;
                negativeIndexes = dataStruct(indShift+j).F_r <=  0;
                x_fit_temp = [dataStruct(indShift+j).t_r, dataStruct(indShift+j).F_r];
                x_fit_temp(negativeIndexes,2) = 1e-25;
                x_fit{j} = x_fit_temp;
                timeVec{j} = dataStruct(indShift+j).t_r;
                maxwellInputs{j} = [dataStruct(indShift+j).t_r(dataStruct(indShift+j).t_r>=0), (dataStruct(indShift+j).h_r(dataStruct(indShift+j).t_r>=0).^1.5)];
                maxwellTimeVec{j} = dataStruct(indShift+j).t_r(dataStruct(indShift+j).t_r>=0);
                x_fit_area{j} = [maxwellInputs{j}(:,1), sqrt(2*r_tip*dataStruct(indShift+j).h_r(dataStruct(indShift+j).t_r>=0) - (dataStruct(indShift+j).h_r(dataStruct(indShift+j).t_r>=0).^2))./r_tip,dataStruct(indShift+j).h_r(dataStruct(indShift+j).t_r>=0)];
            end
        else
            for j = 1:n_datasets
                F_log{j} = log_scale(dataStruct(indShift+j).F_r_smooth,dataStruct(indShift+j).t_r_smooth,dt{j},st{j});
                t_log{j} = log_scale(dataStruct(indShift+j).t_r_smooth,dataStruct(indShift+j).t_r_smooth,dt{j},st{j});
                h_log{j} = log_scale(dataStruct(indShift+j).h_r_smooth,dataStruct(indShift+j).t_r_smooth,dt{j},st{j});
            
                F_schapery_area = dataStruct(indShift+j).F_r_smooth;
                negativeIndexes = dataStruct(indShift+j).F_r_smooth <=  0;
                x_fit_temp = [dataStruct(indShift+j).t_r_smooth, dataStruct(indShift+j).F_r_smooth];
                x_fit_temp(negativeIndexes,2) = 1e-25;
                x_fit{j} = x_fit_temp;
                timeVec{j} = dataStruct(indShift+j).t_r_smooth;
                maxwellInputs{j} = [dataStruct(indShift+j).t_r_smooth(dataStruct(indShift+j).t_r_smooth>=0), (dataStruct(indShift+j).h_r_smooth(dataStruct(indShift+j).t_r_smooth>=0).^1.5)];
                maxwellTimeVec{j} = dataStruct(indShift+j).t_r_smooth(dataStruct(indShift+j).t_r_smooth>=0);
                x_fit_area{j} = [maxwellInputs{j}(:,1), sqrt(2*r_tip*dataStruct(indShift+j).h_r_smooth(dataStruct(indShift+j).t_r_smooth>=0) - (dataStruct(indShift+j).h_r_smooth(dataStruct(indShift+j).t_r_smooth>=0).^2))./r_tip,dataStruct(indShift+j).h_r_smooth(dataStruct(indShift+j).t_r_smooth>=0)];
            end
        end

        % Generate a frequency array in log scale
        de0 = cellfun(@(x)2.0.*pi.*(1./x(end)),timeVec,'UniformOutput',false);
        maxi = cellfun(@(x)2.0.*pi.*(1./mean(diff(x))),timeVec,'UniformOutput',false);
        omega = cellfun(@(x,y)log_tw(x,y),de0,maxi,'UniformOutput',false);

        [dataStruct(1+indShift:endBound).h_norm] = deal(h_norm{:});
        [dataStruct(1+indShift:endBound).F_norm] = deal(F_norm{:});

        [dataStruct(1+indShift:endBound).h_log] = deal(h_log{:});
        h_norm_log = cellfun(@(w,x,y,z)log_scale(w,x,y,z),h_norm,timeVec,dt,st,'UniformOutput',false);
        [dataStruct(1+indShift:endBound).h_norm_log] = deal(h_norm_log{:});
        F_norm_log = cellfun(@(w,x,y,z)log_scale(w,x,y,z),F_norm,timeVec,dt,st,'UniformOutput',false);
        [dataStruct(1+indShift:endBound).F_norm_log] = deal(F_norm_log{:});

        % For Schapery Maxwell
        [dataStruct(1+indShift:endBound).F_schapery_area] = deal(F_schapery_area{:});
        F_log_schapery_area = cellfun(@(w,x,y,z)log_scale(w,x,y,z),F_schapery_area,timeVec,dt,st,'UniformOutput',false);
        [dataStruct(1+indShift:endBound).F_log_schapery_area] = deal(F_log_schapery_area{:});

        % Create the Generalized Voigt model
        [U_func,F_conv,F_conv_wrapper,lb,ub,subref,selector] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);

        % To recreate this data we will need these values later
        
        temp = num2cell(repmat(nu_sample,length(1+indShift:endBound),1));
        [dataStruct(1+indShift:endBound).nu_sample] = deal(temp{:});
        temp = num2cell(repmat(nu_tip,length(1+indShift:endBound),1));
        [dataStruct(1+indShift:endBound).nu_tip] = deal(temp{:});
        temp = num2cell(repmat(E_tip,length(1+indShift:endBound),1));
        [dataStruct(1+indShift:endBound).E_tip] = deal(temp{:});
        temp = num2cell(repmat(r_tip,length(1+indShift:endBound),1));
        [dataStruct(1+indShift:endBound).r_tip] = deal(temp{:});

        % Make Schapery Functions
        [~,E_star,A_conv,h_conv_schapery_area,h_conv_schapery_area_wrapper] = makeSchaperyModel(elasticSetting_maxwell,fluidSetting_maxwell,'n',maxwellTimeVec,dt,st,minTimescale,'maxwell',n_loss_fit,advancedLog,forwardFitTimescale,r_tip,E_tip,nu_tip,nu_sample);

        % Create the Gen. Maxwell Force Convolution Model
        [G_func_maxwell,Q_func_maxwell,h_conv_maxwell,h_conv_maxwell_wrapper,padSize] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_loss_fit,maxwellTimeVec,dt,st);

        % Create the Harmonic Gen. Maxwell Models
        [E_storage,E_loss,lossAngle,harmonic_wrapper,lb_maxwell,ub_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,n_loss_fit,maxwellTimeVec,minTimescale,padSize,forwardFitTimescale);

        % Open Parallel Pool of MATLAB Workers
        if isempty(gcp('nocreate'))
            %parpool(N_workers)
            parpool(8,'IdleTimeout', Inf)
        end

        if timeOperation
            fitSwitchStartTime(iVoigt,1) = toc;
        end

        % Begin fitting
        switch lower(fitMethod)
            case 'iterative'
                linearParams_loop = {};
                maxwellParams_loop = {};
                convParams_loop = {};
                schaperyParams_loop = {};

                if strcmp(elasticSetting,'y')
                    % Get the best elastic parameter for our iterative fit
                    [~,~,~,lb_elastic,ub_elastic,subref_elastic,selector_elastic] = makeGeneralizedVoigtModel('y','n',0,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);

                    x_fit_elastic_mat = cellfun(@(x,y,z)[log_scale(x(:,1),x(:,1),y,z); log_scale(x(:,2),x(:,1),y,z)]',x_fit,dt,st,'UniformOutput',false);

                    F_conv_wrapper_elastic = @(c,input) c(1) .* input(:,2);
                    F_conv_wrapper_elastic_str = '@(c,input) [';
                    % Linear Fit Model Definition
                    for j = 1:n_datasets
                        F_conv_wrapper_elastic_str = horzcat(F_conv_wrapper_elastic_str,sprintf('c(1) .* input{%d}(:,2)',j));
                        if j < n_datasets
                            F_conv_wrapper_elastic_str = horzcat(F_conv_wrapper_elastic_str,';');
                        else
                            F_conv_wrapper_elastic_str = horzcat(F_conv_wrapper_elastic_str,']');
                        end
                    end
                    F_conv_wrapper_elastic_mat = str2func(F_conv_wrapper_elastic_str);

                    y_fit_elastic_mat = cellfun(@transpose,h_norm_log,'UniformOutput',false);

                    if fitElasticFirst
                        beta_dist_elastic = zeros(length(ub_elastic),n_samples);
                        beta0_dist_elastic = zeros(size(beta_dist_elastic));
                        resnorm_dist_elastic = zeros(1,n_samples);

                        if timeOperation
                            elasticParallelStartTime(iVoigt,1) = toc;
                            ticBytes(gcp);
                        end

                        ub_rand_elastic = log10(ub_elastic(1))-1;
                        lb_rand_elastic = log10(lb_elastic(1))+1;
                        beta0_elastic_array = logspace(ub_rand_elastic,lb_rand_elastic,n_samples);

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
                                lsqcurvefit(F_conv_wrapper_elastic_mat,beta0_elastic_array(i),x_fit_elastic_mat,vertcat(y_fit_elastic_mat{:}),lb_elastic,ub_elastic,lsqoptions);

                            beta_dist_elastic(:,i) = betatemp;
                            beta0_dist_elastic(:,i) = beta0_elastic_array(i);

                            if forceResidual
                                for j = 1:n_datasets
                                    yModelIndentation = (1./((8*sqrt(r_tip))/(3*(1-nu_sample))).*F_conv_wrapper_elastic(betatemp,x_fit_elastic_mat{j})).^(2/3);
                                    resnorm_dist_elastic(i) = resnorm_dist_elastic(i)+sum(((yModelIndentation-h_log{j}).^2)./vertcat(movvar(h_log{j},3)))./(length(h_log{j})-length(betatemp)); % Standard Error of Regression (S)
                                end
                            else
                                for j = 1:n_datasets
                                    resnorm_dist_elastic(i) = resnorm_dist_elastic(i)+sum(((F_conv_wrapper_elastic(betatemp,x_fit_elastic_mat{j})-vertcat(y_fit_elastic_mat{j})).^2)./vertcat(movvar(y_fit_elastic_mat{j},3)))./(length(vertcat(y_fit_elastic_mat{j}))-length(betatemp)); % Standard Error of Regression (S)
                                end
                            end
                            
                        end

                        if timeOperation
                            elasticParallelEndTime(iVoigt,1) = toc;
                            elasticParallelBytes(iVoigt,1) = getfield(tocBytes(gcp),{1});

                            fprintf('\nThe Elastic Parameter parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                                (elasticParallelEndTime(iVoigt,1)-elasticParallelStartTime(iVoigt,1))/60, elasticParallelBytes(iVoigt,1)*1e-6);
                        end
                        
                        resnorm_dist_elastic_temp = resnorm_dist_elastic;
                        resnorm_dist_elastic_temp(resnorm_dist_elastic_temp == 0) = NaN;
                        [~,idx] = min(resnorm_dist_elastic_temp,[],'omitnan');
                        bestElasticTerm = beta_dist_elastic(:,idx);
                        
                    else
                        if timeOperation
                            elasticParallelStartTime(iVoigt,1) = toc;
                            ticBytes(gcp);
                        end
                        if timeOperation
                            elasticParallelEndTime(iVoigt,1) = toc;
                            elasticParallelBytes(iVoigt,1) = getfield(tocBytes(gcp),{1});
                        end

                    end
                end

                if timeOperation
                    iterativeFitStartTime(iVoigt,1) = toc;
                end

                % Perform the iterative fitting
                for i = 1:n_terms

                    [~,~,~,~,ub_loop,~,~] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,i,timeVec,dt,st,minTimescale,advancedLog,forwardFitTimescale);
                    [~,~,~,~,padSize_loop] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellTimeVec,dt,st);
                    [~,~,~,~,~,ub_loop_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellTimeVec,minTimescale,padSize_loop,forwardFitTimescale);

                    [tauInds,modulusInds] = getParamIndices(ub_loop,elasticSetting,fluidSetting);
                    [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(ub_loop_maxwell,elasticSetting_maxwell,fluidSetting_maxwell);

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

                    if smoothData == 'n'
                        for j = 1:n_datasets
                            if clipData
                                y_fit_loop{j} = log_scale(h_norm{j}(x_fit{j}(:,1) < max(ub_loop(tauInds))),timeVec{j},dt{j},st{j});
                                y_fit_maxwell_loop{j} = log_scale(F_norm{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec{j},dt{j},st{j});
                                y_fit_area_loop{j} = log_scale(F_schapery_area{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec{j},dt{j},st{j});

                                maxwellInputs_loop{j} = maxwellInputs{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                x_fit_loop{j} = x_fit{j}(x_fit{j}(:,1) < max(ub_loop(tauInds)),:);
                                x_fit_area_loop{j} = x_fit_area{j}(x_fit_area{j}(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                t_log_loop{j} = log_scale(x_fit_loop{j}(:,1),timeVec{j},dt{j},st{j});
                                h_r_loop{j} = dataStruct(indShift+j).h_r(x_fit{j}(:,1) < max(ub_loop(tauInds)));
                                h_norm_loop{j} = h_norm{j}(x_fit{j}(:,1) < max(ub_loop(tauInds)));
                                F_r_loop{j} = dataStruct(indShift+j).F_r(x_fit{j}(:,1) < max(ub_loop(tauInds)));
                                F_norm_loop{j} = F_norm{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                            else
                                y_fit_loop{j} = log_scale(h_norm{j},timeVec{j},dt{j},st{j});
                                y_fit_maxwell_loop{j} = log_scale(F_norm{j},timeVec{j},dt{j},st{j});
                                y_fit_area_loop{j} = log_scale(F_schapery_area{j},timeVec{j},dt{j},st{j});

                                maxwellInputs_loop{j} = maxwellInputs{j};
                                x_fit_loop{j} = x_fit{j};
                                x_fit_area_loop{j} = x_fit_area{j};
                                t_log_loop{j} = log_scale(x_fit_loop{j}(:,1),timeVec{j},dt{j},st{j});
                                h_r_loop{j} = dataStruct(indShift+j).h_r;
                                h_norm_loop{j} = h_norm{j};
                                F_r_loop{j} = dataStruct(indShift+j).F_r;
                                F_norm_loop{j} = F_norm{j};
                            end
                        end
                    else
                        for j = 1:n_datasets
                            if clipData
                                y_fit_loop{j} = log_scale(h_norm{j}(x_fit{j}(:,1) < max(ub_loop(tauInds))),timeVec{j},dt{j},st{j});
                                y_fit_maxwell_loop{j} = log_scale(F_norm{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec{j},dt{j},st{j});
                                y_fit_area_loop{j} = log_scale(F_schapery_area{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds_maxwell))),timeVec{j},dt{j},st{j});

                                maxwellInputs_loop{j} = maxwellInputs{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                x_fit_loop{j} = x_fit{j}(x_fit{j}(:,1) < max(ub_loop(tauInds)),:);
                                x_fit_area_loop{j} = x_fit_area{j}(x_fit_area{j}(:,1) < max(ub_loop_maxwell(tauInds)),:);
                                t_log_loop{j} = log_scale(x_fit_loop{j}(:,1),timeVec{j},dt{j},st{j});
                                h_r_loop{j} = dataStruct(indShift+j).h_r_smooth(x_fit{j}(:,1) < max(ub_loop(tauInds)));
                                h_norm_loop{j} = h_norm{j}(x_fit{j}(:,1) < max(ub_loop(tauInds)));
                                F_r_loop{j} = dataStruct(indShift+j).F_r_smooth(x_fit{j}(:,1) < max(ub_loop(tauInds)));
                                F_norm_loop{j} = F_norm{j}(x_fit{j}(:,1) < max(ub_loop_maxwell(tauInds_maxwell)));
                            else
                                y_fit_loop{j} = log_scale(h_norm{j},timeVec{j},dt{j},st{j});
                                y_fit_maxwell_loop{j} = log_scale(F_norm{j},timeVec{j},dt{j},st{j});
                                y_fit_area_loop{j} = log_scale(F_schapery_area{j},timeVec{j},dt{j},st{j});

                                maxwellInputs_loop{j} = maxwellInputs{j};
                                x_fit_loop{j} = x_fit{j};
                                x_fit_area_loop{j} = x_fit_area{j};
                                t_log_loop{j} = log_scale(x_fit_loop{j}(:,1),timeVec{j},dt{j},st{j});
                                h_r_loop{j} = dataStruct(indShift+j).h_r_smooth;
                                h_norm_loop{j} = h_norm{j};
                                F_r_loop{j} = dataStruct(indShift+j).F_r_smooth;
                                F_norm_loop{j} = F_norm{j};
                            end
                        end
                    end

                    h_norm_log_loop = y_fit_loop;
                    F_norm_log_loop = y_fit_maxwell_loop;
                    h_log_loop = cellfun(@(y,x,d,s)log_scale(y,x(:,1),d,s),h_r_loop,x_fit_loop,dt,st,'UniformOutput',false);
                    F_log_loop = y_fit_area_loop;
                    F_schapery_area_loop = F_r_loop;

                    F_conv_loop = cell(1,n_datasets);
                    F_conv_wrapper_loop = cell(size(F_conv_loop));
                    subref_loop = cell(size(F_conv_loop));
                    h_conv_schapery_area_loop = cell(size(F_conv_loop));
                    h_conv_schapery_area_wrapper_loop = cell(size(F_conv_loop));
                    h_conv_maxwell_loop = cell(size(F_conv_loop));
                    h_conv_maxwell_wrapper_loop = cell(size(F_conv_loop));
                    
                    for j = 1:n_datasets
                        % Make the Generalized Voigt Model for this loop
                        [~,F_conv_loop{j},F_conv_wrapper_loop{j},lb_loop,ub_loop,subref_loop{j},~] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,i,x_fit_loop{j}(:,1),dt{j},st{j},minTimescale,advancedLog,forwardFitTimescale);

                        % Make the Schapery functions for this loop
                        [~,~,~,h_conv_schapery_area_loop{j},h_conv_schapery_area_wrapper_loop{j}] = makeSchaperyModel(elasticSetting_maxwell,fluidSetting_maxwell,'n',maxwellInputs_loop{j}(:,1),dt{j},st{j},minTimescale,'maxwell',i,advancedLog,forwardFitTimescale,r_tip,E_tip,nu_tip,nu_sample);

                        % Create the Gen. Maxwell Force Convolution Model
                        [G_func_maxwell_loop,~,h_conv_maxwell_loop{j},h_conv_maxwell_wrapper_loop{j},padSize_loop] = makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellInputs_loop{j}(:,1),dt{j},st{j});

                        % Create the Harmonic Gen. Maxwell Models
                        [~,~,~,~,lb_loop_maxwell,ub_loop_maxwell] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,i,maxwellInputs_loop{j}(:,1),minTimescale,padSize_loop,forwardFitTimescale); 
                    end
                    
                    beta0_linear = NaN(size(ub_loop));
                    beta0_maxwell = NaN(size(ub_loop_maxwell));
                    beta0_conv = NaN(size(ub_loop));
    %                     beta0_schapery = NaN(size(ub_loop));
                    beta0_schapery = NaN(size(ub_loop_maxwell));

                    ub_randLimit_multigrid = -2*ones(size(ub_loop));
                    ub_randLimit_multigrid(tauInds) = ub_loop(tauInds);

                    lb_randLimit_multigrid = -8*ones(size(lb_loop));
                    lb_randLimit_multigrid(tauInds) = lb_loop(tauInds);

                    ub_randLimit_maxwell_multigrid = 8*ones(size(ub_loop_maxwell));
                    ub_randLimit_maxwell_multigrid(tauInds_maxwell) = ub_loop_maxwell(tauInds_maxwell);

                    lb_randLimit_maxwell_multigrid = 2*ones(size(lb_loop_maxwell));
                    lb_randLimit_maxwell_multigrid(tauInds_maxwell) = lb_loop_maxwell(tauInds_maxwell);

                    ub_randLimit_schapery_multigrid = ub_randLimit_maxwell_multigrid;
                    lb_randLimit_schapery_multigrid = lb_randLimit_maxwell_multigrid;

                    % For the Voigt models
                    if i > 1
                        if fluidSetting == 'y'
                            beta0_linear(1:length(linearParams_loop{i-1})-1) = linearParams_loop{i-1}(1:end-1);
                            beta0_conv(1:length(convParams_loop{i-1})-1) = convParams_loop{i-1}(1:end-1);

                            beta0_linear(end) = 0;
                            beta0_conv(end) = 0;
                        else
                            beta0_linear(1:length(linearParams_loop{i-1})) = linearParams_loop{i-1};
                            if strcmp(elasticSetting,'y') && isnan(beta0_linear(1))
                                beta0_linear(1) = bestElasticTerm;
                            end
                            beta0_conv(1:length(convParams_loop{i-1})) = convParams_loop{i-1};
                            if strcmp(elasticSetting,'y') && isnan(beta0_conv(1))
                                beta0_conv(1) = bestElasticTerm;
                            end
                        end
                    else
                        % Manual override beta0 for random sampling
                        beta0_loop = ((ub_loop-lb_loop).*rand(size(ub_loop)) + lb_loop);
                        beta0_linear = beta0_loop;
                        beta0_conv = beta0_loop;
                    end

                    % For the Maxwell models
                    if i > 1
                        if fluidSetting_maxwell == 'y'
                            beta0_maxwell(1:length(maxwellParams_loop{i-1})-1) = maxwellParams_loop{i-1}(1:end-1);
                            beta0_schapery(1:length(schaperyParams_loop{i-1})-1) = schaperyParams_loop{i-1}(1:end-1);

                            beta0_maxwell(end) = 0;
                            beta0_schapery(end) = 0;
                        else
                            beta0_maxwell(1:length(maxwellParams_loop{i-1})) = maxwellParams_loop{i-1};
                            if strcmp(elasticSetting_maxwell,'y') && isnan(beta0_maxwell(1))
                                beta0_maxwell(1) = 1/bestElasticTerm;
                            end
                            beta0_schapery(1:length(schaperyParams_loop{i-1})) = schaperyParams_loop{i-1};
                            if strcmp(elasticSetting_maxwell,'y') && isnan(beta0_schapery(1))
                                beta0_schapery(1) = 1/bestElasticTerm;
                            end
                        end
                    else
                        % Manual override beta0 for random sampling
                        beta0_loop = ((ub_loop_maxwell-lb_loop_maxwell).*rand(size(ub_loop_maxwell)) + lb_loop_maxwell);
                        beta0_maxwell = beta0_loop;
                        beta0_schapery = beta0_loop;
                    end

                    betaGrid = [];
                    betaGrid_maxwell = [];
                    betaGrid_schapery = [];
                    betaGridIn = {};
                    betaGridIn_maxwell = {};
                    betaGridIn_schapery = {};

                    if bruteForceBeta
                        % Create N-Dimensional Grid of guesses
                        tauStep = 2; % Number of steps in the range of tau starting values (set as low as possible).
                        newInds_conv = isnan(beta0_conv);
                        newInds_maxwell = isnan(beta0_maxwell);
                        newInds_schapery = isnan(beta0_schapery);

                        if strcmp(elasticSetting,'y')
                            elasticInd = 1;
                        else
                            elasticInd = 0;
                        end

                        numStepsPerTerm = floor((n_samples-length(tauInds)*tauStep)^(1./(length(ub_loop(modulusInds))-elasticInd)))/(length(tauInds)*tauStep);

                        for kk = 1:(length(ub_loop))
                            if any(kk == find(newInds_conv,1)) || isempty(find(newInds_conv,1))
                                if (strcmp(elasticSetting,'n') || kk > elasticInd) && fitElasticFirst
                                    if any(kk == tauInds)
                                        loopStep = tauStep;
                                    elseif any(kk == modulusInds)
                                        loopStep = numStepsPerTerm;
                                    end
                                    betaGridIn{kk} = logspace(log10(lb_loop(kk)),log10(ub_loop(kk)),loopStep);
                                else
                                    betaGridIn{kk} = bestElasticTerm;
                                end
                            else
                                betaGridIn{kk} = convParams_loop{i-1};
                            end
                        end

                        for kk = 1:(length(ub_loop_maxwell))
                            if any(kk == find(newInds_maxwell,1)) || isempty(find(newInds_maxwell,1))
                                if (strcmp(elasticSetting_maxwell,'n') || kk > elasticInd) && fitElasticFirst
                                    if any(kk == tauInds_maxwell)
                                        loopStep = tauStep;
                                    elseif any(kk == modulusInds_maxwell)
                                        loopStep = numStepsPerTerm;
                                    end
                                    betaGridIn_maxwell{kk} = logspace(log10(lb_loop_maxwell(kk)),log10(ub_loop_maxwell(kk)),loopStep);
                                    betaGridIn_schapery{kk} = logspace(log10(lb_loop_maxwell(kk)),log10(ub_loop_maxwell(kk)),loopStep);
                                else
                                    betaGridIn_maxwell{kk} = 1./bestElasticTerm;
                                    betaGridIn_schapery{kk} = 1./bestElasticTerm;
                                end
                            else
                                betaGridIn_maxwell{kk} = maxwellParams_loop{i-1};
                                betaGridIn_schapery{kk} = schaperyParams_loop{i-1};
                            end
                        end

                        c = betaGridIn;
                        [c{:}] = ndgrid(c{:});
                        nspecs = length(c);
                        betaGrid = reshape(cat(nspecs+1,c{:}),[],nspecs);

                        c = betaGridIn_maxwell;
                        [c{:}] = ndgrid(c{:});
                        nspecs = length(c);
                        betaGrid_maxwell = reshape(cat(nspecs+1,c{:}),[],nspecs);

                        c = betaGridIn_schapery;
                        [c{:}] = ndgrid(c{:});
                        nspecs = length(c);
                        betaGrid_schapery = reshape(cat(nspecs+1,c{:}),[],nspecs);

                        loopLim = min([size(betaGrid,1), size(betaGrid_maxwell,1), size(betaGrid_schapery,1)]);
                        fprintf('\nBrute Force Grid Size: Iterations Set to %d\n',loopLim);

                    else

                        % Need to create placeholders that look real
                        % because otherwise the parfor pre-processing will
                        % fail.
                        loopLim = n_samples;
                        betaGrid = ones(loopLim,length(ub_loop));
                        betaGrid_maxwell = ones(loopLim,length(ub_loop_maxwell));
                        betaGrid_schapery = ones(loopLim,length(ub_loop_maxwell));

                    end

                    beta_dist = zeros(length(ub_loop),loopLim);
                    beta0_dist = zeros(size(beta_dist));
                    resnorm_dist = zeros(1,loopLim);
                    beta_dist_maxwell = zeros(length(ub_loop_maxwell),loopLim);
                    beta0_dist_maxwell = zeros(size(beta_dist_maxwell));
                    resnorm_dist_maxwell = zeros(1,loopLim);
                    beta_dist_schapery_area = zeros(length(ub_loop_maxwell),loopLim);
                    beta0_dist_schapery_area = zeros(size(beta_dist_schapery_area));
                    resnorm_dist_schapery_area = zeros(1,loopLim);
                    beta_dist_linear = zeros(length(ub_loop),loopLim);
                    beta0_dist_linear = zeros(size(beta_dist_linear));
                    resnorm_dist_linear = zeros(1,loopLim);

                    if timeOperation
                        parallelFitStartTime(iVoigt,1,i) = toc;
                        convTime = Par(loopLim);
                        maxwellTime = Par(loopLim);
                        schaperyTime = Par(loopLim);
                        if i == 1
                           convTimeTotal = [];
                           maxwellTimeTotal = [];
                           schaperyTimeTotal = [];
                        end
                        ticBytes(gcp);
                        if ~exist('userMem','var')
                            %[userMem,systemMem] = memory;
                            userMem = 0;
                            systemMem = 0;
                        end
                        ramused = NaN(1,loopLim);
                        ramav = NaN(1,loopLim);
                        ramtime = NaN(1,loopLim);
                    end

                    errorCalc_str = '@(c,input,y,errorcalc,func) sum([';
                    for j = 1:n_datasets
                        errorCalc_str = horzcat(errorCalc_str,sprintf('errorcalc(c,input{%d},y{%d},func{%d})',j,j,j));
                        if j < n_datasets
                            errorCalc_str = horzcat(errorCalc_str,';');
                        else
                            errorCalc_str = horzcat(errorCalc_str,'])');
                        end
                    end
                    errorCalc_func = str2func(errorCalc_str);
                    
                    errorCalcVoigt = @(c,inputs,y,func) sum(((func(c,inputs)-y).^2)./movvar(y,3))./(length(y)-length(c)); % Standard Error of Regression (S)
                    errorCalcMaxwell = @(c,inputs,y,func) sum(((func(c,inputs)-y).^2)./movvar(y,3))./(length(y)-length(c)); % Standard Error of Regression (S)
                    errorCalcSchapery = @(c,inputs,y,func) sum(((func(c,inputs)-y).^2)./movvar(y,3))./(length(y)-length(c)); % Standard Error of Regression (S)
                    
                    switch lower(optimizationMethod)
                        case 'gd'
                            
                            

                        case 'nls'
                            
                            
                            
                        case {'gen-alg','pattern','particle-swarm','surrogate'}

                            error('Your chosen optimization method is not implemented for Iterative Multi-Load-Level Analysis. Please switch to another optimization method.');

                        case 'nelder-mead'
                            
                            ydata_loop = F_r_loop;
                            ydata_loop_mat = vertcat(ydata_loop{:});
                            ydata_voigt_loop = h_r_loop;
                            ydata_voigt_loop_mat = vertcat(ydata_voigt_loop{:});
                            inputs_loop = cellfun(@(x,y)[x(:,1),y],x_fit_loop,F_r_loop,'UniformOutput',false);
                            inputsForce = cellfun(@(x,y)[x(:,1),y.^1.5],x_fit_loop,h_r_loop,'UniformOutput',false);
                            inputsForce_schapery = cellfun(@(x,y)[x(x(:,1)>0,1),sqrt(2*r_tip*y(x(:,1)>0)-(y(x(:,1)>0).^2))./r_tip],...
                                x_fit_loop,h_r_loop,'UniformOutput',false);
                            timeVec_loop = cellfun(@(x)x(:,1),x_fit_loop,'UniformOutput',false);
                            alphaModel = (8*sqrt(r_tip))/(3*(1-nu_sample));
                            y_fit_GD_loop = h_norm_loop;
                            y_fit_maxwell_GD_loop = F_norm_loop;
                            y_fit_schapery_GD_loop = F_schapery_area_loop;

                            useAnnealing = 0;
                            tossOutOfBounds = 0;
                            plotFit = 0;
                            
                            if plotFit
                                if exist('fitPlot','var')
                                    clearvars fitPlot
                                end
                                if exist('fitPlot','var')
                                    figure(fitPlot);
                                    clf;
                                else
                                    fitPlot = figure('position',[125 25 600 400]);
                                end
                            end

                            parfor j = 1:loopLim
                                warning('off');
                                % Initialize
                                newInds_linear = [];
                                newInds_conv = [];
                                newInds_maxwell = [];
                                newInds_schapery = [];
                                beta0_linear_loop = [];
                                beta0_conv_loop = [];
                                beta0_maxwell_loop = [];
                                beta0_schapery_loop = [];
                                ub_randLimit = [];
                                lb_randLimit = [];
                                ub_randLimit_maxwell = [];
                                lb_randLimit_maxwell = [];
                                ub_randLimit_schapery = [];
                                lb_randLimit_schapery = [];
                                beta0_temp_loop = [];
                                limScale = [];
                                tauInds = [];
                                oldTau = [];
                                modulusInds = [];
                                oldModuli = [];
                                newParams_temp = [];
                                F_conv_loop_weight = [];
                                h_conv_maxwell_loop_weight = [];
                                h_conv_schapery_area_loop_weight = [];
                                beta0_linear_input = [];
                                beta0_linear_input = beta0_linear;
                                beta0_conv_input = [];
                                beta0_conv_input = beta0_conv;
                                beta0_maxwell_input = [];
                                beta0_maxwell_input = beta0_maxwell;
                                beta0_schapery_input = [];
                                beta0_schapery_input = beta0_schapery;

                                % Check for bad indices being fed forward
                                if ( any((ub_loop(~isnan(beta0_linear)) - beta0_linear(~isnan(beta0_linear))) < 0) ) || ( any((beta0_linear(~isnan(beta0_linear)) - lb_loop(~isnan(beta0_linear))) < 0) )
                                    beta0_linear_input(or((ub_loop(~isnan(beta0_linear)) - beta0_linear(~isnan(beta0_linear))) < 0,(beta0_linear(~isnan(beta0_linear)) - lb_loop(~isnan(beta0_linear))) < 0)) = NaN;
                                end
                                if ( any((ub_loop(~isnan(beta0_conv)) - beta0_conv(~isnan(beta0_conv))) < 0) ) || ( any((beta0_conv(~isnan(beta0_conv)) - lb_loop(~isnan(beta0_conv))) < 0) )
                                    beta0_conv_input(or((ub_loop(~isnan(beta0_conv)) - beta0_conv(~isnan(beta0_conv))) < 0,(beta0_conv(~isnan(beta0_conv)) - lb_loop(~isnan(beta0_conv))) < 0)) = NaN;
                                end
                                if ( any((ub_loop_maxwell(~isnan(beta0_maxwell)) - beta0_maxwell(~isnan(beta0_maxwell))) < 0) ) || ( any((beta0_maxwell(~isnan(beta0_maxwell)) - lb_loop_maxwell(~isnan(beta0_maxwell))) < 0) )
                                    beta0_maxwell_input(or((ub_loop_maxwell(~isnan(beta0_maxwell)) - beta0_maxwell(~isnan(beta0_maxwell))) < 0,(beta0_maxwell(~isnan(beta0_maxwell)) - lb_loop_maxwell(~isnan(beta0_maxwell))) < 0)) = NaN;
                                end
                                if ( any((ub_loop_maxwell(~isnan(beta0_schapery)) - beta0_schapery(~isnan(beta0_schapery))) < 0) ) || ( any((beta0_schapery(~isnan(beta0_schapery)) - lb_loop_maxwell(~isnan(beta0_schapery))) < 0) )
                                    beta0_schapery_input(or((ub_loop_maxwell(~isnan(beta0_schapery)) - beta0_schapery(~isnan(beta0_schapery))) < 0,(beta0_schapery(~isnan(beta0_schapery)) - lb_loop_maxwell(~isnan(beta0_schapery))) < 0)) = NaN;
                                end

                                % Find the new indices
                                newInds_linear = isnan(beta0_linear_input);
                                newInds_conv = isnan(beta0_conv_input);
                                newInds_maxwell = isnan(beta0_maxwell_input);
                                newInds_schapery = isnan(beta0_schapery_input);

                                % Control the bounds AND the initial parameters
                                ub_randLimit = ub_randLimit_multigrid;
                                lb_randLimit = lb_randLimit_multigrid;
                                ub_randLimit_maxwell = ub_randLimit_maxwell_multigrid;
                                lb_randLimit_maxwell = lb_randLimit_maxwell_multigrid;
                                ub_randLimit_schapery = ub_randLimit_schapery_multigrid;
                                lb_randLimit_schapery = lb_randLimit_schapery_multigrid;
                                ub_fluidityRandLimit_init = 1e-20;
                                limScale = getfield(logspace(5,0,n_terms),{i})*relaxationFactor;

                                % Include or remove elastic term from bound
                                % loosening that happens below. Use 1 for
                                % inclusion, empty brackets [] for ignoring.
                                elasticInd = 1;

    %                             oldModReductionFactor = limScale*0.5;
                                oldModReductionFactor = 1;
    %                             oldModReductionFactor = 1/relaxationFactor;

                                options = optimset('Display','none',...
                                    'PlotFcns',[],...
                                    'MaxFunEvals',n_maxIterations,...
                                    'MaxIter',n_maxIterations,...
                                    'TolFun',scaling,...
                                    'TolX',0);

                                nelderopts = struct(...
                                    'CoolSched',@(T) (.8*T),...
                                    'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
                                    'InitTemp',1,...
                                    'MaxConsRej',1000,...
                                    'MaxSuccess',20,...
                                    'MaxTries',300,...
                                    'StopTemp',1e-8,...
                                    'StopVal',-Inf,...
                                    'Verbosity',0);

                                % Create random starting point
                                [beta0_conv_loop,tauInds,modulusInds] = makeRandomParams(beta0_conv_input,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,newInds_conv);

                                weightArray = ones(size(beta0_conv_loop));
                                if scaleParams
                                    numTau = length(tauInds);
                                    if strcmp(fluidSetting,'y')
                                       numTau = numTau-1;
                                    end
                                    if i == 1
                                        switch lower(paramScale)
                                            case 'logscale'
                                                F_conv_loop_weight = cellfun(@(x)logscaleParams(x),F_conv_loop,'UniformOutput',false);

                                            otherwise
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(modulusInds) = 10.^((max(ub_randLimit(modulusInds))-min(lb_randLimit(modulusInds)))/2);
                                                F_conv_loop_weight = cellfun(@(x)weightParams(x,weightArray),F_conv_loop,'UniformOutput',false);
                                        end
                                    else
                                        switch lower(paramScale)
                                            case 'equal'
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_conv_input(intersect(find(~newInds_conv),modulusInds))))));
                                                F_conv_loop_weight = cellfun(@(x)weightParams(x,weightArray),F_conv_loop,'UniformOutput',false);

                                            case 'individual'
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(intersect(find(~newInds_conv),modulusInds)) = 10.^(floor(log10(beta0_conv_input(intersect(find(~newInds_conv),modulusInds)))));
                                                weightArray(intersect(find(newInds_conv),modulusInds)) = mean(10.^(floor(log10(beta0_conv_input(intersect(find(~newInds_conv),modulusInds))))));
                                                F_conv_loop_weight = cellfun(@(x)weightParams(x,weightArray),F_conv_loop,'UniformOutput',false);

                                            case 'logscale'
                                                F_conv_loop_weight = cellfun(@(x)logscaleParams(x),F_conv_loop,'UniformOutput',false);
                                        end
                                    end
                                else
                                    F_conv_loop_weight = cellfun(@(x)weightParams(x,weightArray),F_conv_loop,'UniformOutput',false);
                                end

                                if scaleParams
                                    switch lower(paramScale)
                                        case 'logscale'
                                            ub_conv_loop = log10(ub_loop);
                                            lb_conv_loop = log10(lb_loop);
                                            beta0_conv_loop = log10(beta0_conv_loop);
                                        otherwise
                                            ub_conv_loop = (ub_loop)./weightArray;
                                            lb_conv_loop = (lb_loop)./weightArray;
                                            beta0_conv_loop = beta0_conv_loop./weightArray;
                                    end
                                else
                                    ub_conv_loop = (ub_loop);
                                    lb_conv_loop = (lb_loop);
                                end

                                % Gradient Descent
                                if timeOperation
                                    Par.tic;
                                end
                                
%                                 costObj = @(c,inputs,y) sum((F_conv_loop_weight(c,inputs)-y).^2);
                                costObj_str = '@(c,inputs,y,modelfunc) sum([';
                                for k = 1:n_datasets
                                    costObj_str = horzcat(costObj_str,sprintf('((modelfunc{%d}(c,inputs{%d})-y{%d}).^2)',k,k,k));
                                    if k < n_datasets
                                        costObj_str = horzcat(costObj_str,';');
                                    else
                                        costObj_str = horzcat(costObj_str,'])');
                                    end
                                end
                                costObj_func = str2func(costObj_str);
                                costObj_mat = @(c) costObj_func(c,x_fit_loop,y_fit_GD_loop,F_conv_loop_weight);

                                if ~useAnnealing
                                    [newParams_voigt,finalRes,eflag] = fminsearch(@(c) boundObjective(costObj_mat,c,ub_conv_loop,lb_conv_loop), beta0_conv_loop, options);
                                else
                                    [newParams_voigt,finalRes] = annealOpt(@(c) boundObjective(costObj_mat,c,ub_conv_loop,lb_conv_loop), beta0_conv_loop, nelderopts, options);
                                end

                                if ~strcmp(lower(paramScale),'logscale') && scaleParams
                                    newParams_voigt = newParams_voigt.*weightArray;
                                    beta0_conv_loop = beta0_conv_loop.*weightArray;
                                elseif  strcmp(lower(paramScale),'logscale') && scaleParams
                                    newParams_voigt = 10.^(newParams_voigt);
                                    beta0_conv_loop = 10.^(beta0_conv_loop);
                                end

                                if tossOutOfBounds
                                    if ( any((ub_loop - newParams_voigt) < 0) ) || ( any((newParams_voigt - lb_loop) < 0) )
                                        newParams_voigt = NaN(size(newParams_voigt));
                                    end
                                else
                                    if any(newParams_voigt < 0)
                                        newParams_voigt = NaN(size(newParams_voigt));
                                    end
                                end
                                
                                if plotFit
                                    figure(fitPlot)
                                    clf
                                    hold on
                                    cellfun(@(x,y)scatter(x(:,1),y,'bx'),x_fit_loop,y_fit_GD_loop,'UniformOutput',false);
                                    cellfun(@(x,func)plot(x(:,1),func(beta0_conv_loop,x),'b-'),x_fit_loop,F_conv_loop,'UniformOutput',false);
                                    cellfun(@(x,func)plot(x(:,1),func(newParams_voigt,x),'r-'),x_fit_loop,F_conv_loop,'UniformOutput',false);
                                    title('Initial vs. Final Optimized Models [Voigt]')
                                    xlabel('Time [s]')
                                    ylabel('Action Integral Value')
                                    set(gca,'YScale','log')
                                    set(gca,'XScale','log')
                                    legend('Data','Initial','Final')
                                    hold off
                                end

                                beta_dist(:,j) = newParams_voigt;
                                beta0_dist(:,j) = beta0_conv_loop;

                                % Calculate Error Against Data Streams
                                if forceResidual
                                    for k = 1:n_datasets
                                        yModelIndentation = (1./alphaModel.*F_conv_loop{k}(newParams_voigt,inputs_loop{k})).^(2/3);
                                        resnorm_dist(j) = resnorm_dist(j)+sum(((yModelIndentation-ydata_voigt_loop{k}).^2)./movvar(ydata_voigt_loop{k},3))./(length(ydata_voigt_loop{k})-length(newParams_voigt)); % Standard Error of Regression (S)
                                        voigtData{k} = [log_scale(x_fit_loop{k}(:,1),x_fit_loop{k}(:,1),dt{k},st{k})',log_scale(yModelIndentation,x_fit_loop{k}(:,1),dt{k},st{k})'];
                                    end
                                else
                                    resnorm_dist(j) = errorCalc_func(newParams_voigt,x_fit_loop,y_fit_GD_loop,errorCalcVoigt,F_conv_loop); % Standard Error of Regression (S)
%                                     resnorm_dist(j) = finalRes;
                                    voigtData = cellfun(@(x,d,s,func)[log_scale(x(:,1),x(:,1),d,s)',func(newParams_voigt,x)'],x_fit_loop,dt,st,F_conv_wrapper_loop,'UniformOutput',false);
                                end

                                if timeOperation
                                    convTime(j) = Par.toc;
                                end

                                % Create random starting point
                                [beta0_maxwell_loop,tauInds,modulusInds] = makeRandomParams(beta0_maxwell_input,ub_randLimit_maxwell,lb_randLimit_maxwell,elasticSetting_maxwell,fluidSetting_maxwell,newInds_maxwell);

                                weightArray = ones(size(beta0_maxwell_loop));
                                if scaleParams
                                    numTau = length(tauInds);
                                    if strcmp(fluidSetting_maxwell,'y')
                                       numTau = numTau-1;
                                    end
                                    if i == 1
                                        switch lower(paramScale)
                                            case 'logscale'
                                                h_conv_maxwell_loop_weight = cellfun(@(x)logscaleParams(x),h_conv_maxwell_loop,'UniformOutput',false);

                                            otherwise
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(modulusInds) = 10.^((max(ub_randLimit_maxwell(modulusInds))-min(lb_randLimit_maxwell(modulusInds)))/2);
                                                h_conv_maxwell_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_maxwell_loop,'UniformOutput',false);

                                        end
                                    else
                                        switch lower(paramScale)
                                            case 'equal'
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_maxwell_input(~newInds_maxwell)))));
                                                h_conv_maxwell_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_maxwell_loop,'UniformOutput',false);

                                            case 'individual'
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(intersect(find(~newInds_maxwell),modulusInds)) = 10.^(floor(log10(beta0_maxwell_input(intersect(find(~newInds_maxwell),modulusInds)))));
                                                weightArray(intersect(find(newInds_maxwell),modulusInds)) = mean(10.^(floor(log10(beta0_maxwell_input(intersect(find(~newInds_maxwell),modulusInds))))));
                                                h_conv_maxwell_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_maxwell_loop,'UniformOutput',false);

                                            case 'logscale'
                                                h_conv_maxwell_loop_weight = cellfun(@(x)logscaleParams(x),h_conv_maxwell_loop,'UniformOutput',false);

                                        end
                                    end
                                else
                                    h_conv_maxwell_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_maxwell_loop,'UniformOutput',false);
                                end

                                if scaleParams
                                    switch lower(paramScale)
                                        case 'logscale'
                                            ub_maxwell_loop = log10(ub_loop_maxwell);
                                            lb_maxwell_loop = log10(lb_loop_maxwell);
                                            beta0_maxwell_loop = log10(beta0_maxwell_loop);
                                        otherwise
                                            ub_maxwell_loop = (ub_loop_maxwell)./weightArray;
                                            lb_maxwell_loop = (lb_loop_maxwell)./weightArray;
                                            beta0_maxwell_loop = beta0_maxwell_loop./weightArray;
                                    end
                                else
                                    ub_maxwell_loop = (ub_loop_maxwell);
                                    lb_maxwell_loop = (lb_loop_maxwell);
                                end

                                % Gradient Descent
                                if timeOperation
                                    Par.tic;
                                end
                                
%                                 costObj = @(c,inputs,y) sum((h_conv_maxwell_loop_weight(c,inputs)-y).^2);
                                costObj_str = '@(c,inputs,y,modelfunc) sum([';
                                for k = 1:n_datasets
                                    costObj_str = horzcat(costObj_str,sprintf('((modelfunc{%d}(c,inputs{%d})-y{%d}).^2)',k,k,k));
                                    if k < n_datasets
                                        costObj_str = horzcat(costObj_str,';');
                                    else
                                        costObj_str = horzcat(costObj_str,'])');
                                    end
                                end
                                costObj_func = str2func(costObj_str);
                                costObj_mat = @(c) costObj_func(c,maxwellInputs_loop,y_fit_maxwell_GD_loop,h_conv_maxwell_loop_weight);

                                if ~useAnnealing
                                    [newParams_maxwell,finalRes,eflag] = fminsearch(@(c) boundObjective(costObj_mat,c,ub_maxwell_loop,lb_maxwell_loop), beta0_maxwell_loop, options);
                                else
                                    [newParams_maxwell,finalRes] = annealOpt(@(c) boundObjective(costObj_mat,c,ub_maxwell_loop,lb_maxwell_loop), beta0_maxwell_loop, nelderopts, options);
                                end

                                if ~strcmp(lower(paramScale),'logscale') && scaleParams
                                    newParams_maxwell = newParams_maxwell.*weightArray;
                                    beta0_maxwell_loop = beta0_maxwell_loop.*weightArray;
                                elseif strcmp(lower(paramScale),'logscale') && scaleParams
                                    newParams_maxwell = 10.^(newParams_maxwell);
                                    beta0_maxwell_loop = 10.^(beta0_maxwell_loop);
                                end

                                if tossOutOfBounds
                                    if ( any((ub_loop_maxwell - newParams_maxwell) < 0) ) || ( any((newParams_maxwell - lb_loop_maxwell) < 0) )
                                        newParams_maxwell = NaN(size(newParams_maxwell));
                                    end
                                else
                                    if any(newParams_maxwell < 0)
                                        newParams_maxwell = NaN(size(newParams_maxwell));
                                    end
                                end
                                
                                if plotFit
                                    figure(fitPlot)
                                    clf
                                    hold on
                                    cellfun(@(x,y)scatter(x(:,1),y,'bx'),maxwellInputs_loop,y_fit_maxwell_GD_loop,'UniformOutput',false);
                                    cellfun(@(x,func)plot(x(:,1),func(beta0_maxwell_loop,x),'b-'),maxwellInputs_loop,h_conv_maxwell_loop,'UniformOutput',false);
                                    cellfun(@(x,func)plot(x(:,1),func(newParams_maxwell,x),'r-'),maxwellInputs_loop,h_conv_maxwell_loop,'UniformOutput',false);
                                    title('Initial vs. Final Optimized Models [Maxwell]')
                                    xlabel('Time [s]')
                                    ylabel('Action Integral Value')
                                    set(gca,'YScale','log')
                                    set(gca,'XScale','log')
                                    legend('Data','Initial','Final')
                                    hold off
                                end

                                beta_dist_maxwell(:,j) = newParams_maxwell;
                                beta0_dist_maxwell(:,j) = beta0_maxwell_loop;

                                % Calculate Error Against Data Streams
                                if forceResidual
                                    for k = 1:n_datasets
                                        if strcmp(elasticSetting_maxwell,'y')
                                            yModelForce = alphaModel.*(newParams_maxwell(1).*(inputsForce{k}(:,2))+subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce{k}), gradient(inputsForce{k}(:,2),dt{k}),'full'))*dt{k});
                                        else
                                            yModelForce = alphaModel.*(subref_loop(convnfft(G_func_maxwell_loop(newParams_maxwell,inputsForce{k}), gradient(inputsForce{k}(:,2),dt{k}),'full'))*dt{k});
                                        end
                                        resnorm_dist_maxwell(j) = resnorm_dist_maxwell(j)+sum(((yModelForce-ydata_loop{k}).^2)./movvar(ydata_loop{k},3))./(length(ydata_loop{k})-length(newParams_maxwell)); % Standard Error of Regression (S)
                                        maxwellData = [log_scale(maxwellInputs_loop{k}(:,1),maxwellInputs_loop{k}(:,1),dt{k},st{k})',log_scale(yModelForce,maxwellInputs_loop{k}(:,1),dt{k},st{k})'];
                                    end
                                else
                                    resnorm_dist_maxwell(j) = errorCalc_func(newParams_maxwell,maxwellInputs_loop,y_fit_maxwell_GD_loop,errorCalcMaxwell,h_conv_maxwell_loop); % Standard Error of Regression (S)
%                                     resnorm_dist_maxwell(j) = finalRes;
                                    maxwellData = cellfun(@(x,d,s,func)[log_scale(x(:,1),x(:,1),d,s)',func(newParams_maxwell,x)'],maxwellInputs_loop,dt,st,h_conv_maxwell_wrapper_loop,'UniformOutput',false);
                                end

                                if timeOperation
                                    maxwellTime(j) = Par.toc;
                                end

                                % Create random starting point
                                [beta0_schapery_loop,tauInds,modulusInds] = makeRandomParams(beta0_schapery_input,ub_randLimit_schapery,lb_randLimit_schapery,elasticSetting_maxwell,fluidSetting_maxwell,newInds_schapery);

                                weightArray = ones(size(beta0_schapery_loop));
                                if scaleParams
                                    numTau = length(tauInds);
                                    if strcmp(fluidSetting_maxwell,'y')
                                       numTau = numTau-1;
                                    end
                                    if i == 1
                                        switch lower(paramScale)
                                            case 'logscale'
                                                h_conv_schapery_area_loop_weight = cellfun(@(x)logscaleParams(x),h_conv_schapery_area_loop,'UniformOutput',false);

                                            otherwise
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(modulusInds) = 10.^((max(ub_randLimit_schapery(modulusInds))-min(lb_randLimit_schapery(modulusInds)))/2);
                                                h_conv_schapery_area_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_schapery_area_loop,'UniformOutput',false);

                                        end
                                    else
                                        switch lower(paramScale)
                                            case 'equal'
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(modulusInds) = 10.^(floor(log10(mean(beta0_schapery_input(~newInds_schapery)))));
                                                h_conv_schapery_area_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_schapery_area_loop,'UniformOutput',false);

                                            case 'individual'
                                                weightArray(tauInds(1:numTau)) = logspace(minTimescale,minTimescale*(10^numTau),numTau);
                                                weightArray(intersect(find(~newInds_schapery),modulusInds)) = 10.^(floor(log10(beta0_schapery_input(intersect(find(~newInds_schapery),modulusInds)))));
                                                weightArray(intersect(find(newInds_schapery),modulusInds)) = mean(10.^(floor(log10(beta0_schapery_input(intersect(find(~newInds_schapery),modulusInds))))));
                                                h_conv_schapery_area_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_schapery_area_loop,'UniformOutput',false);

                                            case 'logscale'
                                                h_conv_schapery_area_loop_weight = cellfun(@(x)logscaleParams(x),h_conv_schapery_area_loop,'UniformOutput',false);

                                        end
                                    end
                                else
                                    h_conv_schapery_area_loop_weight = cellfun(@(x)weightParams(x,weightArray),h_conv_schapery_area_loop,'UniformOutput',false);
                                end

                                if scaleParams
                                    switch lower(paramScale)
                                        case 'logscale'
                                            ub_schapery_loop = log10(ub_loop_maxwell);
                                            lb_schapery_loop = log10(lb_loop_maxwell);
                                            beta0_schapery_loop = log10(beta0_schapery_loop);
                                        otherwise
                                            ub_schapery_loop = (ub_loop_maxwell)./weightArray;
                                            lb_schapery_loop = (lb_loop_maxwell)./weightArray;
                                            beta0_schapery_loop = beta0_schapery_loop./weightArray;
                                    end
                                else
                                    ub_schapery_loop = (ub_loop_maxwell);
                                    lb_schapery_loop = (lb_loop_maxwell);
                                end

                                % Gradient Descent
                                if timeOperation
                                    Par.tic;
                                end

%                                 costObj = @(c,inputs,y) sum((h_conv_schapery_area_loop_weight(c,inputs)-y).^2);
                                costObj_str = '@(c,inputs,y,modelfunc) sum([';
                                for k = 1:n_datasets
                                    costObj_str = horzcat(costObj_str,sprintf('((modelfunc{%d}(c,inputs{%d})-y{%d}).^2)',k,k,k));
                                    if k < n_datasets
                                        costObj_str = horzcat(costObj_str,';');
                                    else
                                        costObj_str = horzcat(costObj_str,'])');
                                    end
                                end
                                costObj_func = str2func(costObj_str);
                                costObj_mat = @(c) costObj_func(c,x_fit_area_loop,y_fit_schapery_GD_loop,h_conv_schapery_area_loop_weight);

                                if ~useAnnealing
                                    [newParams_schapery,finalRes,eflag] = fminsearch(@(c) boundObjective(costObj_mat,c,ub_schapery_loop,lb_schapery_loop), beta0_schapery_loop, options);
                                else
                                    [newParams_schapery,finalRes] = annealOpt(@(c) boundObjective(costObj_mat,c,ub_schapery_loop,lb_schapery_loop), beta0_schapery_loop, nelderopts, options);
                                end

                                if ~strcmp(lower(paramScale),'logscale') && scaleParams
                                    newParams_schapery = newParams_schapery.*weightArray;
                                    beta0_schapery_loop = beta0_schapery_loop.*weightArray;
                                elseif strcmp(lower(paramScale),'logscale') && scaleParams
                                    newParams_schapery = 10.^(newParams_schapery);
                                    beta0_schapery_loop = 10.^(beta0_schapery_loop);
                                end

                                if tossOutOfBounds
                                    if ( any((ub_loop_maxwell - newParams_schapery) < 0) ) || ( any((newParams_schapery - lb_loop_maxwell) < 0) )
                                        newParams_schapery = NaN(size(newParams_schapery));
                                    end
                                else
                                    if any(newParams_schapery < 0)
                                        newParams_schapery = NaN(size(newParams_schapery));
                                    end
                                end
                                
                                if plotFit
                                    figure(fitPlot)
                                    clf
                                    hold on
                                    cellfun(@(x,y)scatter(x(:,1),y,'bx'),x_fit_area_loop,y_fit_schapery_GD_loop,'UniformOutput',false);
                                    cellfun(@(x,func)plot(x(:,1),func(beta0_schapery_loop,x),'b-'),x_fit_area_loop,h_conv_schapery_area_loop,'UniformOutput',false);
                                    cellfun(@(x,func)plot(x(:,1),func(newParams_schapery,x),'r-'),x_fit_area_loop,h_conv_schapery_area_loop,'UniformOutput',false);
                                    title('Initial vs. Final Optimized Models [Schapery]')
                                    xlabel('Time [s]')
                                    ylabel('Action Integral Value')
                                    set(gca,'YScale','log')
                                    set(gca,'XScale','log')
                                    legend('Data','Initial','Final')
                                    hold off
                                end

                                beta_dist_schapery_area(:,j) = newParams_schapery;
                                beta0_dist_schapery_area(:,j) = beta0_schapery_loop;

                                % Calculate Error Against Data Streams
                                if forceResidual
                                    for k = 1:n_datasets
                                        yModelForce = h_conv_schapery_area_loop{k}(newParams_schapery,inputsForce_schapery{k});
                                        resnorm_dist_schapery_area(j) = resnorm_dist_schapery_area(j)+sum(((yModelForce-ydata_loop{k}).^2)./movvar(ydata_loop{k},3))./(length(ydata_loop{k})-length(newParams_schapery)); % Standard Error of Regression (S)
                                        schaperyData{k} = [log_scale(x_fit_area_loop{k}(:,1),x_fit_area_loop{k}(:,1),dt{k},st{k})',log_scale(yModelForce,x_fit_area_loop{k}(:,1),dt{k},st{k})'];
                                    end
                                else
                                    resnorm_dist_schapery_area(j) = errorCalc_func(newParams_schapery,x_fit_area_loop,y_fit_schapery_GD_loop,errorCalcSchapery,h_conv_schapery_area_loop); % Standard Error of Regression (S)
%                                     resnorm_dist_schapery_area(j) = finalRes;
                                    schaperyData = cellfun(@(x,d,s,func)[log_scale(x(:,1),x(:,1),d,s)',func(newParams_schapery,x)'],x_fit_area_loop,dt,st,h_conv_schapery_area_wrapper_loop,'UniformOutput',false);
                                end

                                if timeOperation
                                    schaperyTime(j) = Par.toc;
                                end

                                if timeOperation
                                    % Track memory
                                    %ramused(j) = userMem.MemUsedMATLAB;
                                    ramused(j) = 0;
                                    %ramav(j) = systemMem.PhysicalMemory.Available;
                                    ramav(j) = 0;
                                    ramtime(j) = now;
                                end

                                % Placeholder for the linear model since we
                                % don't care about it.
                                beta_dist_linear(:,j) = ones(size(beta0_conv_loop));
                                beta0_dist_linear(:,j) = ones(size(beta0_conv_loop));
                                resnorm_dist_linear(j) = 1;

                            end % End Parfor
                            
                        otherwise
                            
                            error('That optimization method does not exist.');

                    end
                    
                    if timeOperation
                        parallelFitEndTime(iVoigt,1,i) = toc;
                        parallelFitParallelBytes(iVoigt,1,i) = getfield(tocBytes(gcp),{1});

                        stop(convTime);
                        stop(maxwellTime);
                        stop(schaperyTime);

                        convTimeTotal = [convTimeTotal convTime];
                        maxwellTimeTotal = [maxwellTimeTotal maxwellTime];
                        schaperyTimeTotal = [schaperyTimeTotal schaperyTime];

                        fprintf('\nThe parallel loop took %5.2f minutes, and on average %4.4g Mb were sent to each worker.\n',...
                            (parallelFitEndTime(iVoigt,1,i)-parallelFitStartTime(iVoigt,1,i))/60, parallelFitParallelBytes(iVoigt,1,i)*1e-6);
                    end

                    resnorm_dist_linear_temp = resnorm_dist_linear;
                    resnorm_dist_linear_temp(resnorm_dist_linear_temp == 0) = NaN;
                    [~,idx_linear] = min(resnorm_dist_linear_temp,[],'omitnan');
                    linearParams_loop{i} = beta_dist_linear(:,idx_linear);

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

                    resnorm_dist_schapery_area_temp = resnorm_dist_schapery_area;
                    resnorm_dist_schapery_area_temp(resnorm_dist_schapery_area_temp == 0) = NaN;
                    resnorm_dist_schapery_area_temp( any(abs(log10(beta_dist_schapery_area(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                    resnorm_dist_schapery_area_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))) - log10(beta_dist_schapery_area(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                    [~,idx_schapery] = min(resnorm_dist_schapery_area_temp,[],'omitnan');
                    schaperyParams_loop{i} = beta_dist_schapery_area(:,idx_schapery);
                    
                end

                if timeOperation
                    iterativeFitEndTime(iVoigt,1) = toc;
                    parallelFitBreakdown{iVoigt,1} = {convTimeTotal,maxwellTimeTotal,schaperyTimeTotal};
                    parallelFitMemoryMonitor{iVoigt,1} = vertcat(ramused,ramav,ramtime);
                end

                temp = repmat({linearParams_loop},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).linearParams_loop] = deal(temp{:});
                temp = repmat({convParams_loop},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).convParams_loop] = deal(temp{:});
                temp = repmat({schaperyParams_loop},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).schaperyParams_loop] = deal(temp{:});
                temp = repmat({maxwellParams_loop},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).maxwellParams_loop] = deal(temp{:});

                n_ranks = 10; % Number of places to save the parameter sets of

                resnorm_dist_linear_temp = resnorm_dist_linear;
                resnorm_dist_linear_temp(resnorm_dist_linear_temp == 0) = NaN;
                [~,idx_linear] = min(resnorm_dist_linear_temp,[],'omitnan');
                
                temp = repmat({beta_dist_linear(:,idx_linear)},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).linearParams] = deal(temp{:});

                resnorm_dist_temp = resnorm_dist;
                resnorm_dist_temp(resnorm_dist_temp == 0) = NaN;
                resnorm_dist_temp( any(abs(log10(beta_dist(modulusInds,:)) - log10(repmat(lb_loop(modulusInds)',1,size(beta_dist(modulusInds,:),2))))<tossThresh,1) ) = NaN;
                resnorm_dist_temp( any(abs(log10(repmat(ub_loop(modulusInds)',1,size(beta_dist(modulusInds,:),2))) - log10(beta_dist(modulusInds,:)))<tossThresh,1) ) = NaN;
                [~,idx] = min(resnorm_dist_temp,[],'omitnan');
                temp = repmat({beta_dist(:,idx)},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).convParams] = deal(temp{:});

                [~,I] = sort(resnorm_dist_temp,'MissingPlacement','last');
                temp = repmat({beta_dist(:,I(1:n_ranks))},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).convParams_alternatives] = deal(temp{:});

                resnorm_dist_maxwell_temp = resnorm_dist_maxwell;
                resnorm_dist_maxwell_temp(resnorm_dist_maxwell_temp == 0) = NaN;
                resnorm_dist_maxwell_temp( any(abs(log10(beta_dist_maxwell(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_maxwell(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                resnorm_dist_maxwell_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_maxwell(modulusInds_maxwell,:),2))) - log10(beta_dist_maxwell(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                [~,idx_maxwell] = min(resnorm_dist_maxwell_temp);
                temp = repmat({beta_dist_maxwell(:,idx_maxwell)},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).maxwellParams] = deal(temp{:});

                [~,I_maxwell] = sort(resnorm_dist_maxwell_temp,'MissingPlacement','last');
                temp = repmat({beta_dist_maxwell(:,I_maxwell(1:n_ranks))},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).maxwellParams_alternatives] = deal(temp{:});

                resnorm_dist_schapery_area_temp = resnorm_dist_schapery_area;
                resnorm_dist_schapery_area_temp(resnorm_dist_schapery_area_temp == 0) = NaN;
                resnorm_dist_schapery_area_temp( any(abs(log10(beta_dist_schapery_area(modulusInds_maxwell,:)) - log10(repmat(lb_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))))<tossThresh,1) ) = NaN;
                resnorm_dist_schapery_area_temp( any(abs(log10(repmat(ub_loop_maxwell(modulusInds_maxwell)',1,size(beta_dist_schapery_area(modulusInds_maxwell,:),2))) - log10(beta_dist_schapery_area(modulusInds_maxwell,:)))<tossThresh,1) ) = NaN;
                [~,idx_schapery] = min(resnorm_dist_schapery_area_temp,[],'omitnan');
                temp = repmat({beta_dist_schapery_area(:,idx_schapery)},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).schaperyParams] = deal(temp{:});

                [~,I_schapery] = sort(resnorm_dist_schapery_area_temp,'MissingPlacement','last');
                temp = repmat({beta_dist_schapery_area(:,I_schapery(1:n_ranks))},length(1+indShift:endBound),1);
                [dataStruct(1+indShift:endBound).schaperyParams_alternatives] = deal(temp{:});

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
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                else
                    fprintf(fid,'No elastic response used.\r\n');

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,idx_maxwell),iii,beta_dist_maxwell(iii*2,idx_maxwell));
                    end
                    if fluidSetting_maxwell == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,idx_maxwell));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                end

                fprintf(fid,'\r\n\r\nSchapery Generalized Maxwell\r\n');
                if elasticSetting_maxwell == 'y'
                    fprintf(fid,'E_e = %3.6g [Pa]\r\n',beta_dist_schapery_area(1,idx_schapery));

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,idx_schapery),iii,beta_dist_schapery_area(iii*2+1,idx_schapery));
                    end
                    if fluidSetting_maxwell == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
                    else
                        fprintf(fid,'No fluidity response used.\r\n');
                    end

                else
                    fprintf(fid,'No elastic response used.\r\n');

                    for iii = 1:n_terms
                        fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,idx_schapery),iii,beta_dist_schapery_area(iii*2,idx_schapery));
                    end
                    if fluidSetting_maxwell == 'y'
                        fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,idx_schapery));
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
                        fprintf(fid,'E_g = %3.6g [Pa]\r\n',beta_dist_maxwell(1,I_maxwell(alti)));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2+1,I_maxwell(alti)));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_maxwell(iii*2-1,I_maxwell(alti)),iii,beta_dist_maxwell(iii*2,I_maxwell(alti)));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_maxwell(end,I_maxwell(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                    % Print Schapery Maxwell
                    fprintf(fid,'\r\nSchapery Generalized Maxwell\n');
                    if elasticSetting_maxwell == 'y'
                        fprintf(fid,'E_g = %3.6g [Pa]\r\n',beta_dist_schapery_area(1,I_schapery(alti)));

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2+1,I_schapery(alti)));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    else
                        fprintf(fid,'No elastic response used.\r\n');

                        for iii = 1:n_terms
                            fprintf(fid,'E_%d = %3.6g [Pa], tau_%d = %3.6g [s]\r\n',iii,beta_dist_schapery_area(iii*2-1,I_schapery(alti)),iii,beta_dist_schapery_area(iii*2,I_schapery(alti)));
                        end
                        if fluidSetting_maxwell == 'y'
                            fprintf(fid,'phi_f = %3.6g [1/(Pa*s)]\r\n',beta_dist_schapery_area(end,I_schapery(alti)));
                        else
                            fprintf(fid,'No fluidity response used.\r\n');
                        end

                    end

                end
                        
            case 'open'

                error('Open Search is not implemented for Multi-Load-Level Analysis. Please switch to Iterative.');

        end
        
        if timeOperation
            fitSwitchEndTime(iVoigt,1) = toc;
        end

        if strcmp(avAnalysis,'i')
            save([pathname sprintf('/%s-dataStruct-%d_terms-%s.mat',savePrepend,n_terms,saveLabel)],'dataStruct')
            fprintf('\nSaved the dataStruct to an m-file for %d terms on File No. %d\n',n_terms,i_loop);
        end

        if saveGDdata && strcmp(optimizationMethod,'gd')
            % Set plotting threshold
            plotThreshold = Inf; % Standard Error (S)
            bigFontSize = 28;
            mediumFontSize = 16;
            markersize = 100;
            nPlotsPerRow = 3;
            nRow = ceil(size(beta_dist,1)/nPlotsPerRow);

            for ij = 1:3

                X_beta = [];
                Cost_beta = [];
                alphaBeta = [];
                modelLabel = [];

                switch ij
                    case 1
                        X_beta = (beta_dist(:,resnorm_dist<plotThreshold)');
                        Cost_beta = repmat(resnorm_dist(resnorm_dist<plotThreshold)',[1 size(beta_dist,1)]).*(1./alphaModel);
                        alphaBeta = (min(resnorm_dist(resnorm_dist<plotThreshold))./resnorm_dist(resnorm_dist<plotThreshold))';
                        modelLabel = 'GKV';

                        % Make labels
                        [tauInds,modulusInds] = getParamIndices(X_beta,elasticSetting,fluidSetting);
                        symLabels = cell(size(X_beta,1),1);
                        if strcmp(fluidSetting,'y')
                            symLabels(tauInds) = horzcat(sprintfc('\\tau_{%d}',(1:(length(tauInds)-1))),sprintf('\\phi_{f}'));
                        else
                            symLabels(tauInds) = sprintfc('\\tau_{%d}',(1:(length(tauInds))));
                        end
                        if strcmp(elasticSetting,'y')
                            symLabels(modulusInds) = horzcat(sprintf('J_{g}'),sprintfc('J_{%d}',(1:(length(modulusInds)-1))));
                        else
                            symLabels(modulusInds) = {sprintfc('J_{%d}',(1:(length(modulusInds))))};
                        end

                        yLabelString = 'Standard Error, S [$m^{3/2}$]';

                    case 2
                        X_beta = (beta_dist_maxwell(:,resnorm_dist_maxwell<plotThreshold)');
                        Cost_beta = repmat(resnorm_dist_maxwell',[1 size(beta_dist_maxwell,1)]).*alphaModel;
                        alphaBeta = (min(resnorm_dist_maxwell(resnorm_dist_maxwell<plotThreshold))./resnorm_dist_maxwell(resnorm_dist_maxwell<plotThreshold))';
                        modelLabel = 'GM';

                        % Make labels
                        [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(X_beta,elasticSetting_maxwell,fluidSetting_maxwell);
                        symLabels = cell(size(X_beta,1),1);
                        if strcmp(fluidSetting_maxwell,'y')
                            symLabels(tauInds) = horzcat(sprintfc('\\tau_{%d}',(1:(length(tauInds)-1))),sprintf('\\phi_{f}'));
                        else
                            symLabels(tauInds) = sprintfc('\\tau_{%d}',(1:(length(tauInds))));
                        end
                        if strcmp(elasticSetting_maxwell,'y')
                            symLabels(modulusInds) = horzcat(sprintf('E_{e}'),sprintfc('E_{%d}',(1:(length(modulusInds)-1))));
                        else
                            symLabels(modulusInds) = {sprintfc('E_{%d}',(1:(length(modulusInds))))};
                        end

                        yLabelString = 'Standard Error, S [$N$]';

                    case 3
                        X_beta = (beta_dist_schapery_area(:,resnorm_dist_schapery_area<plotThreshold)');
                        Cost_beta = repmat(resnorm_dist_schapery_area',[1 size(beta_dist_schapery_area,1)]).*alphaModel;
                        alphaBeta = (min(resnorm_dist_schapery_area(resnorm_dist_schapery_area<plotThreshold))./resnorm_dist_schapery_area(resnorm_dist_schapery_area<plotThreshold))';
                        modelLabel = 'SGM';

                        % Make labels
                        [tauInds_maxwell,modulusInds_maxwell] = getParamIndices(X_beta,elasticSetting_maxwell,fluidSetting_maxwell);
                        symLabels = cell(size(X_beta,1),1);
                        if strcmp(fluidSetting_maxwell,'y')
                            symLabels(tauInds) = horzcat(sprintfc('\\tau_{%d}',(1:(length(tauInds)-1))),sprintf('\\phi_{f}'));
                        else
                            symLabels(tauInds) = sprintfc('\\tau_{%d}',(1:(length(tauInds))));
                        end
                        if strcmp(elasticSetting_maxwell,'y')
                            symLabels(modulusInds) = horzcat(sprintf('E_{e}'),sprintfc('E_{%d}',(1:(length(modulusInds)-1))));
                        else
                            symLabels(modulusInds) = {sprintfc('E_{%d}',(1:(length(modulusInds))))};
                        end

                        yLabelString = 'Standard Error, S [$N$]';

                end

                % Begin visualizations of best parameter sets
                if exist('paramSpacePlot','var')
                    clearvars paramSpacePlot
                end
                if exist('paramSpacePlot','var')
                    figure(paramSpacePlot);
                    clf;
                else
                    paramSpacePlot = figure('position',[100 100 1200 500*nRow]);
                end

                for ii = 1:size(beta_dist,1)
                    if ii <= nPlotsPerRow
                        subtightplot (nRow, nPlotsPerRow, ii, [0.1 0.05], [0.2 0.1], [0.075 0.05]);
                    elseif ii > nPlotsPerRow && ii < nPlotsPerRow*(nRow)-(nPlotsPerRow-1)
                        subtightplot (nRow, nPlotsPerRow, ii, [0.1 0.05], [0.2 0.075], [0.075 0.05]);
                    else
                        subtightplot (nRow, nPlotsPerRow, ii, [0.1 0.05], [0.2 0.075], [0.075 0.05]);
                    end
                    hold on
                    for iii = 1:length(alphaBeta)
                        try
%                                 scatter(Cost_beta(iii,ii),X_beta(iii,ii),markersize,'filled','MarkerFaceAlpha',alphaBeta(iii))
                            scatter(X_beta(iii,ii),Cost_beta(iii,ii),markersize,'filled','MarkerFaceAlpha',alphaBeta(iii))
                        catch ERROR
                        end
                    end
                    title(sprintf('$%s$',symLabels{ii}),'interpreter','latex', 'FontSize', bigFontSize)
                    if size(beta_dist,1)-ii < nPlotsPerRow
%                             xlabel({'' '' 'Standard Error, S [$m^{3/2}$]'}, 'interpreter','latex', 'FontSize', mediumFontSize)
                        xlabel({'' '' 'Parameter Value [Various]'}, 'interpreter', 'latex', 'FontSize', mediumFontSize)
                    end
                    if ~mod(ii-1,nPlotsPerRow)
%                                 ylabel('Parameter Value [Various]', 'interpreter', 'latex', 'FontSize', mediumFontSize)
                        ylabel(yLabelString, 'interpreter','latex', 'FontSize', mediumFontSize)
                    end
                    grid on
                    grid minor
                    set(gca,'XScale','log')
                    set(gca,'YScale','log')
                    hold off
                end

                saveas(paramSpacePlot, [pathname sprintf('/ParameterSpaceGD-nElements_%d-%s-%s',n_terms,modelLabel,saveLabel)], 'jpg')
                savefig(paramSpacePlot, [pathname sprintf('/ParameterSpaceGD-nElements_%d-%s-%s',n_terms,modelLabel,saveLabel)],'compact')

                close(paramSpacePlot);

            end

            if strcmp(avAnalysis,'a')
                save([pathname sprintf('/%s-FittingInfo-%d_terms-average.mat',savePrepend,n_terms)],...
                    'beta0_dist','beta0_dist_maxwell','beta0_dist_schapery_area',...
                    'beta_dist','beta_dist_maxwell','beta_dist_schapery_area',...
                    'resnorm_dist','resnorm_dist_maxwell','resnorm_dist_schapery_area')
                fprintf('\nSaved the Fitting Data to an m-file for %d terms and all velocities combined.\n',n_terms);
            end

            if exist('paramSpacePlot','var')
                clearvars paramSpacePlot
            end

        end

        fprintf(fid,'\r\n\r\n=================\r\n\r\n');
        
    end
    
    if timeOperation
        loopEndTimeStore(iVoigt) = toc;
    end
    
    fclose('all');
    
    if strcmp(avAnalysis,'a')
        save([pathname sprintf('/%s-dataStruct-%d_terms-average.mat',savePrepend,n_terms)],'dataStruct')
        fprintf('\nSaved the dataStruct to an m-file for %d terms and all velocities.\n',n_terms);
    end
            
end % End iVoigt

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