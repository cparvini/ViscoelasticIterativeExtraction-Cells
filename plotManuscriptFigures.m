clear all
close all
clc

addpath(genpath(pwd))

% ======================== %
% Viscoelastic Parameter
% Extraction, Iterative
% Methodology Manuscript
%
% Figure Generation
%
% Created by cparvini
% ======================== %
% This file will create the 
% figures presented in the 
% manuscript: "Viscoelastic
% Parameterization of Human
% Skin Cells Characterize 
% Material Behaviour at 
% Multiple Timescales", by
% Parvini, Cartagena-Rivera,
% and Solares, originally
% submitted to Nature Comm.
% Biology in 2021.
%
% To execute this script, 
% one must have run the 
% analysis of the AFM force
% curves in the "data" 
% directory using the 
% "VicsoParameterExtraction.m"
% script.
%
% The analyzed results from the
% manuscript can be provided
% under reasonable request, 
% but are excluded from
% the Github repository due
% to size constraints. In 
% total, the analyzed data 
% amounts to ~1.5 GB.
%
% ======================== %

%% Load MAT File and Parse
getfolderString = 'Please select a folder to add to the comparison';
list = {'GKV','GM'};
modelSelected = list{2};

% Force Plot Choices
dataLabel = {'Melanocytes','Melanoma','Fibroblasts'};
hertzPlotLabel = {'Melanocytes','Melanomas','Fibroblasts'};

% Here, select which cell and which "load level" (i.e. curve number) to use
% for the example plots.
cellSelected = [12 2 1];
loadLevelSelected = [1 4 1];

pathname = cell(size(dataLabel));
for i = 1:length(dataLabel)
    if i > 1
        startdir = pathname{i-1};
        temp = strfind(startdir,{'\'});
        startdir = startdir(1:temp(end));
    else
        startdir = pwd;
    end
    pathname{i} = uigetdir(startdir,...
            sprintf('Please select the parent directory containing the %s cell results',dataLabel{i}));
end

% Include Special External Datasets
% Use a list of strings to include each paper's dataset. Current Options:

% Alcaraz2003 - "Microrheology of Human Lung Epithelial Cells Measured by
% Atomic Force Microscopy" - PLR

plotExternal = {'Alcaraz2003'};

if ~exist('outputpathname','var')
    startdir = pathname{end};
    temp = strfind(startdir,{'\'});
    startdir = startdir(1:temp(end));
    outputpathname = uigetdir(startdir,...
    'Select the Output Folder For Your Results');
end
% saveLabel = input('What would you like the labels to be for your output plots?: ','s');
saveLabel = 'IterativePaper3';

%% Settings and Placeholders
errorCalc = @(ydata,ymodel) sum((ydata-ymodel).^2); % SSE
SError = @(ydata,ymodel) sqrt( sum( ((ydata-ymodel).^2) ./ (numel(ydata) - 2) ) ); % Standard Error

bigFontSize = 40;
mediumFontSize = 34;
smallFontSize = 24;
linewidthSize = 5;
markerSize = 25;
xScaling = 'log';
yScaling = 'log';
boxSetting = 'on'; % 'on' or 'off'
CI = 1; % 1 or 0
CI_method = 'studentt'; % Choose 'observed' or 'studentt'
alphaCI = 0.05;
fillAlpha = 0.2;
harmonicConfig = 'separate'; % use 'single' 'stack' or 'separate'
splitHarmonicRange = 0;
yyLossAngle = false;
if yyLossAngle
    
    leftPlotNums = [1];
    rightPlotNums = [2 3];
    
    for i = 1:length(pathname)
        if ismember(leftPlotNums,i)
            fprintf('The %s Cells will be plotted on the LEFT side of the Loss Angle Plot.\n\n',dataLabel{i});
        else
            fprintf('The %s Cells will be plotted on the RIGHT side of the Loss Angle Plot.\n\n',dataLabel{i});
        end
    end
    
end
useLegend = 0;
limMargin = 0.25;
plotMin = 1;    % Minimum plot limit, in Pa, for the Moduli
plotHertz = true;
colorSpecList = {'r' 'b' 'g' 'm' 'k' 'c'}';
lineSpecList = {'-' '--' ':' '-.'}';
markerSpecList = {'' 'd' 'o' 's' '*' '+' 'v' '>' 'h'}';

close all

if exist('dataFig','var')
    clearvars dataFig
end
if exist('forceFig','var')
    clearvars forceFig
end
if exist('harmonicFig','var')
    clearvars harmonicFig
end
if exist('lossAngleFig','var')
    clearvars lossAngleFig
end

if useLegend
    dataFig = figure('position',[25 25 1600 1000]);
    forceFig = figure('position',[75 25 1600 800]);
    indFig = figure('position',[75 25 1600 800]);
    if ~strcmp(harmonicConfig,'separate')
        harmonicFig = figure('position',[100 25 1600 1000]);
    else
        harmonicFig_storage = figure('position',[100 25 1600 1000]);
        harmonicFig_loss = figure('position',[100 50 1600 1000]);
    end
    if ~isempty(plotExternal)
        lossAngleFig = figure('position',[125 25 1600 1000]);
    else
        lossAngleFig = figure('position',[125 25 1000 1000]);
    end
    hertzFig = figure('position',[150 25 1000 1000]);
    hertzForceFig = figure('position',[175 25 1000 1000]);
else
    dataFig = figure('position',[25 25 1000 1000]);
    forceFig = figure('position',[75 25 1000 500]);
    indFig = figure('position',[75 25 1000 500]);
    if ~strcmp(harmonicConfig,'separate')
        harmonicFig = figure('position',[100 25 1000 1000]);
    else
        harmonicFig_storage = figure('position',[100 25 1000 1000]);
        harmonicFig_loss = figure('position',[100 50 1000 1000]);
    end
    if ~isempty(plotExternal)
        lossAngleFig = figure('position',[125 25 1600 1000]);
    else
        lossAngleFig = figure('position',[125 25 1000 1000]);
    end
    hertzFig = figure('position',[150 25 1000 1000]);
    hertzForceFig = figure('position',[175 25 1000 1000]);
end

legendString = dataLabel;

if exist('dataStruct','var')
    clearvars dataStruct
end

figure(forceFig)
clf
hold on
grid on
set(gca,'xscale',xScaling)
set(gca,'yscale',yScaling)
set(gca,'Box',boxSetting)
set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize)
xlabel('Time [$s$]', 'Interpreter', 'latex', 'FontSize', mediumFontSize)
ylabel('Force [$N$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
hold off

figure(indFig)
clf
hold on
grid on
set(gca,'xscale',xScaling)
set(gca,'yscale',yScaling)
set(gca,'Box',boxSetting)
set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize)
xlabel('Time [$s$]', 'Interpreter', 'latex', 'FontSize', mediumFontSize)
ylabel('Indentation [$m$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
hold off

if ~strcmp(harmonicConfig,'separate')
    figure(harmonicFig)
    clf
else
    figure(harmonicFig_storage)
    clf
    figure(harmonicFig_loss)
    clf
end

figure(lossAngleFig)
clf
hold on
lossAngOrig = findobj(get(gcf, 'Children'), '-depth', 1, 'type', 'axes');
colorOrder = get(lossAngOrig, 'ColorOrder');
grid on
set(gca,'Box',boxSetting)
set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize)
set(gca,'xscale',xScaling)
xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
if yyLossAngle
    yyaxis left
    set(gca,'yscale','linear')
    ylabel('$\delta$($\omega$) [$deg$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
    set(gca,'YColor',[0 0 0])
    yyaxis right
    set(gca,'yscale','linear')
    ylabel('$\delta$($\omega$) [$deg$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
    set(gca,'YColor',[0 0 0])
else
    set(gca,'yscale','linear')
    ylabel('$\delta$($\omega$) [$deg$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
end
set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
hold off

figure(hertzFig)
clf
hold on
set(gca,'Box',boxSetting)
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize)
ylabel('Young''s Modulus [$kPa$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
hold off

figure(hertzForceFig)
clf
hold on
grid on
set(gca,'xscale',xScaling)
set(gca,'yscale',yScaling)
set(gca,'Box',boxSetting)
set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize)
xlabel('Time [$s$]', 'Interpreter', 'latex', 'FontSize', mediumFontSize)
ylabel('Force [$N$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
hold off

%% Loop through the directories
for i_dir = 1:length(pathname)
    
    % Initialize the data figure for this cell type
    figure(dataFig)
    clf
    hold on
    grid on
    set(gca,'xscale',xScaling)
    set(gca,'yscale',yScaling)
    set(gca,'Box',boxSetting)
    set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
    set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize)
    xlabel('Time [$s$]', 'Interpreter', 'latex', 'FontSize', mediumFontSize)
    switch modelSelected
        case 'GKV'
            ylabel('$\int_{0}^{t} U(t-\xi) F(\xi) d\xi$ [$m^2$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
        case 'GM'
            ylabel('$\int_{0}^{t} Q(t-\xi) h^{3/2}(\xi) d\xi$ [$N/\sqrt(m)$]', 'FontSize', mediumFontSize, 'Interpreter', 'Latex');
    end
    hold off
    dataPlot_legendString = {};
    
    % Check to see if there are subdirectories
    dirContents = dir(pathname{i_dir});
    subFolders = dirContents([dirContents.isdir]);
    subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];
    subFolders = natsortfiles(subFolders);
    
    % If the user provides a main directory with many subdirectories containing
    % data, we should loop through all directories and analyze each in turn.
    if ~isempty(subFolders)
        Folders = cell(1,length(subFolders));
        Folders = cellfun(@(root,sub)[root '\' sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
    else
        Folders = pathname(i_dir);
    end
    
    % Placeholder for the optimal number of terms for each velocity in each
    % directory
    optimalTermNumbers = cell(size(Folders));
    optimalParams = cell(size(Folders));
    optimalHertz = cell(size(Folders));
    optimalOmegas = cell(size(Folders));
    
    % Find the optimal performing models for each cell
    for i = 1:length(Folders)
        % Analyze each directory in the cell folders
        % Find the settings file in this folder
        settingsFile = dir([Folders{i} filesep '*settingsStruct*.mat']);
        load([settingsFile.folder filesep settingsFile.name]);
        
        avAnalysis = settingsStruct.avAnalysis;
        elasticSetting = settingsStruct.elasticSetting;
        fluidSetting = settingsStruct.fluidSetting;
        elasticSetting_maxwell = settingsStruct.elasticSetting_maxwell;
        fluidSetting_maxwell = settingsStruct.fluidSetting_maxwell;
        minTimescale = settingsStruct.minTimescale;
        forwardFitTimescale = settingsStruct.forwardFitTimescale;

        % First, we must determine the "best fit" number of terms for this
        % directory. There may be multiple velocities, in which case they
        % should be treated separately for determining the "best number of
        % terms".
        Files = dir([Folders{i} '\*dataStruct*.mat']);
        errorMat = Inf(length(Files),1);
        hertzErrorMat = Inf(length(Files),1);
        paramList = cell(size(errorMat));
        hertzModulus = cell(size(hertzErrorMat));
        
        for j = 1:length(Files)
            
            % Load the file
            clearvars dataStruct
            load([Files(j).folder '\' Files(j).name]);
            
            % Find the number of terms
            temp = strsplit(Files(j).name,{'-','.'});
            termind = 3;
            if strcmpi(temp{1},'NIH')
                termind = termind+1;
            end
            temp = strsplit(temp{termind},{'_'});
            n_terms = str2num(temp{1});
            
            % Calculate the offsets
            indShift = 0;
            for k = 1:size(dataStruct,2)
                if isempty(dataStruct(k).t_average)
                    indShift = indShift + 1;
                end
            end
            
            % Find the number of averages. Make sure we have enough room in
            % our error matrix.
            n_averages = size(dataStruct,2)-indShift;
            if strcmp(avAnalysis,'i')
                loopMax = indShift;
                indShift = 0;
                
                if size(errorMat,2) < loopMax
                    errorMat = horzcat(errorMat,repmat(Inf(length(Files),1),1,loopMax-size(errorMat,2)));
                    hertzErrorMat = horzcat(hertzErrorMat,repmat(Inf(length(Files),1),1,loopMax-size(hertzErrorMat,2)));
                end
                
            else
                loopMax = n_averages;
                
                if size(errorMat,2) < n_averages
                    errorMat = horzcat(errorMat,repmat(Inf(length(Files),1),1,n_averages-size(errorMat,2)));
                    hertzErrorMat = horzcat(hertzErrorMat,repmat(Inf(length(Files),1),1,n_averages-size(hertzErrorMat,2)));
                end
                
            end
            
            r_tip = settingsStruct.r_tip;
            nu_sample = settingsStruct.nu_sample;

            alpha = (8*sqrt(r_tip))/(3*(1-nu_sample));
            
            F_hertz = @(E,hdata) (4*sqrt(r_tip)/(3*(1-nu_sample^2))).*E.*(hdata.^1.5);
                        
            % Loop through the files in this directory
            for k = 1:loopMax
                switch lower(modelSelected)
                    case 'gm'
                        if settingsStruct.smoothData == 'n'
                            inputs = [dataStruct(indShift+k).t_r,...
                                dataStruct(indShift+k).h_r.^1.5];
                        else
                            inputs = [dataStruct(indShift+k).t_r_smooth,...
                                dataStruct(indShift+k).h_r_smooth.^1.5];
                        end
                        timeVec = inputs(:,1);
                        dt = dataStruct(indShift+k).dt;
                        st = timeVec(end);
                        [~,~,h_conv,~,~] = ...
                            makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,...
                            n_terms,timeVec,dt,st);
                        modelParams = dataStruct(indShift+k).maxwellParams;
                        xData = timeVec;
                        yData = dataStruct(indShift+k).F_norm;
                        xModel = timeVec;
                        yModel = h_conv(modelParams,inputs);
                        fitError = errorCalc(yData,yModel);
                    case 'gkv'
                        if settingsStruct.smoothData == 'n'
                            inputs = [dataStruct(indShift+k).t_r,...
                                dataStruct(indShift+k).F_r];
                        else
                            inputs = [dataStruct(indShift+k).t_r_smooth,...
                                dataStruct(indShift+k).F_r_smooth];
                        end
                        timeVec = inputs(:,1);
                        dt = dataStruct(indShift+k).dt;
                        st = timeVec(end);
                        [~,F_conv,~,~,~,~,~] = ...
                            makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,...
                            timeVec,dt,st,minTimescale,forwardFitTimescale);
                        modelParams = dataStruct(indShift+k).convParams;
                        xData = timeVec;
                        yData = dataStruct(indShift+k).h_norm;
                        xModel = timeVec;
                        yModel = F_conv(modelParams,inputs);
                        fitError = errorCalc(yData,yModel);
                end
                errorMat(j,k) = fitError;
                paramList{j,k} = modelParams;
                
                % Grab the Hertzian Modulus
                hertzModulus{j,k} = dataStruct(indShift+k).hertzianModulus;
                if settingsStruct.smoothData == 'n'
                    tdata = dataStruct(indShift+k).t_r;
                    hdata = dataStruct(indShift+k).h_r;
                    yData = dataStruct(indShift+k).F_r;
                else
                    tdata = dataStruct(indShift+k).t_r_smooth;
                    hdata = dataStruct(indShift+k).h_r_smooth;
                    yData = dataStruct(indShift+k).F_r_smooth;
                end
                yModel = F_hertz(hertzModulus{j,k},hdata);
                hertzErrorMat(j,k) = errorCalc(yData,yModel);
                
            end
            
        end
        
        errorMargin = 0.01;
        [~,bestidx] = min(errorMat,[],1);
        [~,bestidxHertz] = min(hertzErrorMat,[],1);
        bestErrors = Inf(size(bestidx));
        for j = 1:length(bestidx)
            bestErrors(j) = errorMat(bestidx(j),j);
        end
        closeErrors = abs((errorMat-repmat(bestErrors,size(errorMat,1),1))./(repmat(bestErrors,size(errorMat,1),1))) < errorMargin;
        for j = 1:length(bestidx)
            if isnan(bestErrors(j))
                continue;
            end
            bestidx(j) = find(closeErrors(:,j),1,'first');
            bestErrors(j) = errorMat(bestidx(j),j);
        end
        
        optimalTermNumbers{i} = bestidx;
        for j = 1:length(bestidx)
            
            if isnan(bestErrors(j))
                continue;
            end
            
            optimalParams{i} = horzcat(optimalParams{i},paramList(bestidx(j),j));
            optimalHertz{i} = horzcat(optimalHertz{i},hertzModulus(bestidxHertz(j),j));
            
            if plotHertz && i == cellSelected(i_dir) && j == loadLevelSelected(i_dir)
                if settingsStruct.smoothData == 'n'
                    tdata = dataStruct(indShift+j).t_r;
                    hdata = dataStruct(indShift+j).h_r;
                    yData = dataStruct(indShift+j).F_r;
                else
                    tdata = dataStruct(indShift+j).t_r_smooth;
                    hdata = dataStruct(indShift+j).h_r_smooth;
                    yData = dataStruct(indShift+j).F_r_smooth;
                end
                yModel = F_hertz(hertzModulus{bestidxHertz(j),j},hdata);

                xModelHertz = tdata;
                yModelHertz = yModel;
                xDataHertz = tdata;
                yDataHertz = yData;
            end
            
            % Grab the info we will need below from the optimal set
            % Load the file
            clearvars dataStruct
            load([Files(bestidx(j)).folder '\' Files(bestidx(j)).name]);

            % Find the number of terms
            temp = strsplit(Files(bestidx(j)).name,{'-','.'});
            termind = 3;
            if strcmpi(temp{1},'NIH')
                termind = termind+1;
            end
            temp = strsplit(temp{termind},{'_'});
            n_terms = str2num(temp{1});

            % Calculate the offsets
            indShift = 0;
            for k = 1:size(dataStruct,2)
                if isempty(dataStruct(k).t_average)
                    indShift = indShift + 1;
                end
            end
            
            n_averages = size(dataStruct,2)-indShift;
            if strcmp(avAnalysis,'i')
                indShift = 0;
            end
            
            switch lower(modelSelected)
                case 'gm'
                    if settingsStruct.smoothData == 'n'
                        inputs = [dataStruct(indShift+j).t_r,...
                            dataStruct(indShift+j).h_r.^1.5];
                    else
                        inputs = [dataStruct(indShift+j).t_r_smooth,...
                            dataStruct(indShift+j).h_r_smooth.^1.5];
                    end
                    timeVec = inputs(:,1);
                    dt = dataStruct(indShift+j).dt;
                    st = timeVec(end);
                    [~,~,h_conv,~,~] = ...
                        makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,...
                        n_terms,timeVec,dt,st);
                case 'gkv'
                    if settingsStruct.smoothData == 'n'
                        inputs = [dataStruct(indShift+j).t_r,...
                            dataStruct(indShift+j).F_r];
                    else
                        inputs = [dataStruct(indShift+j).t_r_smooth,...
                            dataStruct(indShift+j).F_r_smooth];
                    end
                    timeVec = inputs(:,1);
                    dt = dataStruct(indShift+j).dt;
                    st = timeVec(end);
                    [~,F_conv,~,~,~,~,~] = ...
                        makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,...
                        timeVec,dt,st,minTimescale,forwardFitTimescale);
            end
            de0 = 2.0.*pi.*(1./st);
            maxi = 2.0.*pi.*(1./dt);
            omega = log_tw(de0,maxi);
            optimalOmegas{i} = horzcat(optimalOmegas{i},{omega});
        end
        
    end
    
    % Data Plot
    % Figure 4, Subfigures c, d, and e
    % First, load the data we are choosing to plot
    % Load the file
    Files = dir([Folders{cellSelected(i_dir)} '\*dataStruct*.mat']);
        
    for j = 1:length(Files)
        styleSpec = [colorSpecList{j}];
        
        linewidthSize_loop = linewidthSize*getfield(linspace(2.5,1,length(Files)),{j});
%         linewidthSize_loop = linewidthSize;
        
        % Load the file
        clearvars dataStruct
        load([Files(j).folder '\' Files(j).name]);
            
        % Find the number of terms
        temp = strsplit(Files(j).name,{'-','.'});
        termind = 3;
        if strcmpi(temp{1},'NIH')
            termind = termind+1;
        end
        temp = strsplit(temp{termind},{'_'});
        n_terms = str2num(temp{1});
            
        % Calculate the offsets
        indShift = 0;
        for k = 1:size(dataStruct,2)
            if isempty(dataStruct(k).t_average)
                indShift = indShift + 1;
            end
        end
        if strcmp(avAnalysis,'i')
            indShift = 0;
        end

        % Loop through the files in this directory
        switch lower(modelSelected)
            case 'gm'
                if settingsStruct.smoothData == 'n'
                    inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r,...
                        dataStruct(indShift+loadLevelSelected(i_dir)).h_r.^1.5];
                else
                    inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth,...
                        dataStruct(indShift+loadLevelSelected(i_dir)).h_r_smooth.^1.5];
                end
                timeVec = inputs(:,1);
                dt = dataStruct(indShift+loadLevelSelected(i_dir)).dt;
                st = timeVec(end);
                [~,~,h_conv,~,~] = ...
                    makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,...
                    n_terms,timeVec,dt,st);
                modelParams = dataStruct(indShift+loadLevelSelected(i_dir)).maxwellParams;
                xData = timeVec;
                yData = dataStruct(indShift+loadLevelSelected(i_dir)).F_norm;
                xModel = timeVec;
                yModel = h_conv(modelParams,inputs);
                xModelForce = timeVec;
                yModelForce = alpha.*h_conv(modelParams,inputs);
                xDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).t_r;
                yDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).F_r;
                
            case 'gkv'
                if settingsStruct.smoothData == 'n'
                    inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r,...
                        dataStruct(indShift+loadLevelSelected(i_dir)).F_r];
                else
                    inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth,...
                        dataStruct(indShift+loadLevelSelected(i_dir)).F_r_smooth];
                end
                timeVec = inputs(:,1);
                dt = dataStruct(indShift+loadLevelSelected(i_dir)).dt;
                st = timeVec(end);
                [~,F_conv,~,~,~,~,~] = ...
                    makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,...
                    timeVec,dt,st,minTimescale,forwardFitTimescale);
                modelParams = dataStruct(indShift+loadLevelSelected(i_dir)).convParams;
                xData = timeVec;
                yData = dataStruct(indShift+loadLevelSelected(i_dir)).h_norm;
                xModel = timeVec;
                yModel = F_conv(modelParams,inputs);
                
                % Convert Parameters for calculating force
                if settingsStruct.smoothData == 'n'
                    inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r,...
                        dataStruct(indShift+loadLevelSelected(i_dir)).h_r.^1.5];
                else
                    inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth,...
                        dataStruct(indShift+loadLevelSelected(i_dir)).h_r_smooth.^1.5];
                end
                [~,~,h_conv,~,padSize] = ...
                    makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,...
                    n_terms,timeVec,dt,st);
                newparams = convertParamsCollocation(modelParams,'gm',elasticSetting,fluidSetting,elasticSetting_maxwell,fluidSetting_maxwell,n_terms,padSize,timeVec,dt,st,minTimescale,forwardFitTimescale);
                xModelForce = timeVec;
                yModelForce = alpha.*h_conv(modelParams,inputs);
                xDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).t_r;
                yDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).F_r;
        end

        % Plot the arm comparison
        figure(dataFig)
        hold on
        if j == 1
            scatter(xData,yData,markerSize*10,'bx','Linewidth',linewidthSize/2)
            dataPlot_legendString{j} = 'Cell Data';
        end
        plot(xModel,yModel,styleSpec,'LineWidth',linewidthSize_loop);
        hold off
        dataPlot_legendString{j+1} = sprintf('%d Term(s)', n_terms);
        
        fprintf('Data Figure, %s, %d Term(s), Standard Error: %4.4g\n',dataLabel{i_dir},n_terms,SError(yData,yModel))
        
    end
    
    % Save Data Plot
    figure(dataFig)
    limSet = findobj(gca, '-property', 'ydata');
    limSet = get(limSet, 'YData');
    limSet = [limSet{:}];
    set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
    if strcmp(get(findobj(gca, '-property', 'xscale'),'XScale'),'log')
        limSet = findobj(gca, '-property', 'xdata');
        limSet = get(limSet, 'XData');
        limSet = [limSet{:}];
        set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
    end
    saveas(dataFig, [outputpathname '\' sprintf('%s-DataFitPlot-%s',saveLabel,dataLabel{i_dir})], 'jpg')
    saveas(dataFig, [outputpathname '\' sprintf('%s-DataFitPlot-%s',saveLabel,dataLabel{i_dir})], 'fig')
    
    if useLegend && (i_dir == 1)
        legendFig = figure;
        hold on
        % Plot whatever you like
        x = 1:10;
        y = NaN;
        for i = 1:length(Files)
            plot(x, y .* x, 'DisplayName', sprintf('%d Term(s)',i),'color',colorSpecList{i},'lineWidth',linewidthSize)
        end
        % Initial values to capture the entire legend
        % Should fit most modern screens
        set(gcf,'Position',[0,0,1024,1024]);
        % Call the legend to your choice, I used a horizontal legend here
        legend_handle = legend('Orientation','vertical','FontSize',mediumFontSize);
        % Set the figure Position using the normalized legend Position vector
        % as a multiplier to the figure's current position in pixels
        % This sets the figure to have the same size as the legend
        set(gcf,'Position',(get(legend_handle,'Position')...
            .*[0, 0, 1, 1].*get(gcf,'Position')));
        % The legend is still offset so set its normalized position vector to
        % fill the figure
        set(legend_handle,'Position',[0,0,1,1]);
        % Put the figure back in the middle screen area
        set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);

        saveas(legendFig, [outputpathname '\' sprintf('%s-DataFitPlotLegend',saveLabel)], 'jpg')
        saveas(legendFig, [outputpathname '\' sprintf('%s-DataFitPlotLegend',saveLabel)], 'fig')
        close(legendFig)
        clearvars legendFig
    end
    
    % Harmonic Plots
    % Figures 5 and 6
    % Analyze the frequencies available from each cell, and choose only the
    % range that ALL datasets can reasonably provide information on.
    w_min = 0;
    w_max = Inf;
    for i = 1:length(optimalOmegas)
        for j = 1:length(optimalOmegas{i})
            if min(optimalOmegas{i}{j}) > w_min
                w_min = min(optimalOmegas{i}{j});
            end
            if max(optimalOmegas{i}{j}) < w_max
                w_max = max(optimalOmegas{i}{j});
            end
        end
    end
    
    % Create our "universal" omega frequency array for plotting
    omega = logspace(log10(w_min),log10(w_max),10*abs(log10(w_max)-log10(w_min)));
    if ~isrow(omega)
        omega = omega';
    end
    xLossMod = omega;
    xStorageMod = omega;
    xLossAngle = omega;
    
    % NOW we can loop through all of our parameter sets and calculate the
    % harmonics from every cell
    lossModAll = [];
    storageModAll = [];
    lossAngleAll = [];
    for i = 1:length(optimalParams)
        for j = 1:length(optimalParams{i})
            switch lower(modelSelected)
                case 'gm'
                    n_terms = (length(optimalParams{i}{j})-strcmp(fluidSetting_maxwell,'y')-strcmp(elasticSetting_maxwell,'y'))/2;
                    [~,~,~,~,padSize] = ...
                            makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,...
                            n_terms,timeVec,dt,st);
                    [E_storage,E_loss,lossAngle,~,~,~] = ...
                            makeHarmonicGeneralizedMaxwellModel(fluidSetting_maxwell,...
                            elasticSetting_maxwell,n_terms,timeVec,minTimescale,...
                            padSize,forwardFitTimescale);
                    % Get the harmonics
                    storageMod = E_storage(omega,optimalParams{i}{j});
                    lossMod = E_loss(omega,optimalParams{i}{j});
                    lossAngle = lossAngle(omega,optimalParams{i}{j});
                case 'gkv'
                    Jabs = sqrt(J_storage_advanced(omega,optimalParams{i}{j},elasticSetting,fluidSetting).^2 + ...
                        J_loss_advanced(omega,optimalParams{i}{j},elasticSetting,fluidSetting).^2);
                    storageMod = J_storage_advanced(omega,optimalParams{i}{j},elasticSetting,fluidSetting)./(Jabs.^2);
                    lossMod = J_loss_advanced(omega,optimalParams{i}{j},elasticSetting,fluidSetting)./(Jabs.^2);
                    lossAngle = atand(lossMod./storageMod);
            end
            
            % Stack the harmonics
            storageModAll = vertcat(storageModAll,storageMod);
            lossModAll = vertcat(lossModAll,lossMod);
            lossAngleAll = vertcat(lossAngleAll,lossAngle);
            
        end
    end
    
    % We also need to loop through all of our hertzian moduli and collect
    % them for later.
    if ~exist('hertzList','var')
        hertzList = [];
        hertzLabel = {};
    end
    
    for i = 1:length(optimalHertz)
        temp = cell2mat(optimalHertz{i});
        if isrow(temp) temp = temp'; end
        hertzList = vertcat(hertzList,temp);
        hertzLabel = vertcat(hertzLabel,repmat(hertzPlotLabel(i_dir),length(temp),1));
    end
    
    % Here we can create a plot showing ALL the harmonics from ALL the
    % curves, if we want!
    
    % Now, we create the harmonic confidence intervals
    if CI
        
        CI_storage = NaN(2,size(storageModAll,2));
        CI_loss = NaN(2,size(lossModAll,2));
        CI_angle = NaN(2,size(lossAngleAll,2));
        
        switch lower(CI_method)
            case 'observed'
                % Find the harmonic bounds represented by our datasets.
                CI_storage(1,:) = max(storageModAll,[],1);
                CI_storage(2,:) = min(storageModAll,[],1);
                CI_loss(1,:) = max(lossModAll,[],1);
                CI_loss(2,:) = min(lossModAll,[],1);
                CI_angle(1,:) = max(lossAngleAll,[],1);
                CI_angle(2,:) = min(lossAngleAll,[],1);
                
            case 'studentt'
                % Use the Student T table to calculate the range of
                % expected values for each timestep.
                dof = size(storageModAll,1);
                CI_storage(1,:) = mean(storageModAll,1) + tinv(1-alphaCI/2,dof).*std(storageModAll,0,1)./sqrt(size(storageModAll,1));
                CI_storage(2,:) = mean(storageModAll,1) + tinv(alphaCI/2,dof).*std(storageModAll,0,1)./sqrt(size(storageModAll,1));
                
                dof = size(lossModAll,1);
                CI_loss(1,:) = mean(lossModAll,1) + tinv(1-alphaCI/2,dof).*std(lossModAll,0,1)./sqrt(size(lossModAll,1));
                CI_loss(2,:) = mean(lossModAll,1) + tinv(alphaCI/2,dof).*std(lossModAll,0,1)./sqrt(size(lossModAll,1));
                
                dof = size(lossAngleAll,1);
                CI_angle(1,:) = mean(lossAngleAll,1) + tinv(1-alphaCI/2,dof).*std(lossAngleAll,0,1)./sqrt(size(lossAngleAll,1));
                CI_angle(2,:) = mean(lossAngleAll,1) + tinv(alphaCI/2,dof).*std(lossAngleAll,0,1)./sqrt(size(lossAngleAll,1));
                
            otherwise
                error('You have selected an unknown Confidence Interval method. Please check your selection.');
        end
        
        % Create our "averages"
        yStorageMod = mean(storageModAll,1);
        yLossMod = mean(lossModAll,1);
        yLossAngle = mean(lossAngleAll,1);
        
        CI_storage = sort(CI_storage,1,'descend');
        CI_loss = sort(CI_loss,1,'descend');
        CI_angle = sort(CI_angle,1,'descend');
        
        % Remove negative numbers
        CI_storage(CI_storage < 0) = plotMin;
        CI_loss(CI_loss < 0) = plotMin;
        CI_angle(CI_angle < 0) = plotMin;
        
    end
    
    % Harmonic Plot
    switch harmonicConfig
        case 'single'
            figure(harmonicFig)
            yyaxis left
            hold on
            plot(xStorageMod,yStorageMod,'Linewidth',linewidthSize_loop)
            if CI
                set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                h1a = plot(xStorageMod,CI_storage(1,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                h1b = plot(xStorageMod,CI_storage(2,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
            end
            grid on
            set(gca,'xscale',xScaling)
            set(gca,'yscale',yScaling)
            set(gca,'ylim',yLims)
            set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop)
            xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
            ylabel('Storage Modulus [$Pa$]', 'FontSize', bigFontSize, 'Interpreter', 'Latex');

            hold off

            yyaxis right
            hold on
            plot(xLossMod,yLossMod,'Linewidth',linewidthSize_loop)
            if CI
                set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                h1a = plot(xLossMod,CI_loss(1,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                h1b = plot(xLossMod,CI_loss(2,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
            end
            grid on
            set(gca,'xscale',xScaling)
            set(gca,'yscale',yScaling)
            set(gca,'ylim',yLims)
            set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop)
            xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
            ylabel('Loss Modulus [$Pa$]', 'FontSize', bigFontSize, 'Interpreter', 'Latex');
            hold off

        case 'stack'
            figure(harmonicFig)
            if ~splitHarmonicRange
                subtightplot (2, 1, 1, [0.05 0.05], [0.125 0.05], [0.125 0.05]);
                hold on
                grid on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
        %         set(findall(gcf,'-property','FontSize'),'FontSize',smallFontSize)
                set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop/2)
        %         xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
                ylabel('Storage Modulus [$Pa$]', 'FontSize', smallFontSize, 'Interpreter', 'Latex');
                plot(xStorageMod,yStorageMod,'Linewidth',linewidthSize_loop)
                if CI
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1a = plot(xStorageMod,CI_storage(1,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1b = plot(xStorageMod,CI_storage(2,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                end
                set(gca,'xticklabel',[])
                if i > 1
                    limSet = findobj(gca, '-property', 'ydata');
                    limSet = get(limSet, 'YData');
                    if iscell(limSet)
                        limSet = [limSet{:}];
                    end
                    set(gca,'ylim',[min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)])
                end
                hold off

                subtightplot (2, 1, 2, [0.05 0.05], [0.125 0.05], [0.125 0.05]);
                hold on
                grid on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                set(findall(gcf,'-property','FontSize'),'FontSize',smallFontSize)
                set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop/2)
                xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', smallFontSize)
                ylabel('Loss Modulus [$Pa$]', 'FontSize', smallFontSize, 'Interpreter', 'Latex');
                plot(xLossMod,yLossMod,'Linewidth',linewidthSize_loop)
                if CI
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1a = plot(xLossMod,CI_loss(1,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1b = plot(xLossMod,CI_loss(2,:),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                end
                if i > 1
                    limSet = findobj(gca, '-property', 'ydata');
                    limSet = get(limSet, 'YData');
                    if iscell(limSet)
                        limSet = [limSet{:}];
                    end
                    set(gca,'ylim',[min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)])
                end
                hold off
            else
                subtightplot (2, 2, 1, [0.05 0.05], [0.125 0.05], [0.125 0.05]);
                hold on
                grid on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
    %             set(findall(gcf,'-property','FontSize'),'FontSize',smallFontSize)
                set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop/2)
    %             xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
                ylabel('Storage Modulus [$Pa$]', 'FontSize', smallFontSize, 'Interpreter', 'Latex');
                plot(xStorageMod(1:round(length(xStorageMod)/2)),yStorageMod(1:round(length(xStorageMod)/2)),'Linewidth',linewidthSize_loop)
                if CI
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1a = plot(xStorageMod(1:round(length(xStorageMod)/2)),CI_storage(1,(1:round(length(xStorageMod)/2))),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1b = plot(xStorageMod(1:round(length(xStorageMod)/2)),CI_storage(2,(1:round(length(xStorageMod)/2))),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                end
                set(gca,'xlim',[xStorageMod(1) xStorageMod(round(length(xStorageMod)/2))])
                set(gca,'xticklabel',[])
    %             ylim auto
                hold off
                if i > 1
                    limSet = findobj(gca, '-property', 'ydata');
                    limSet = get(limSet, 'YData');
                    if iscell(limSet)
                        limSet = [limSet{:}];
                    end
                    set(gca,'ylim',[min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)])
                end
                
                subtightplot (2, 2, 2, [0.05 0.05], [0.125 0.05], [0.15 0.025]);
                hold on
                grid on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
    %             set(findall(gcf,'-property','FontSize'),'FontSize',smallFontSize)
                set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop/2)
    %             xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
    %             ylabel('Storage Modulus [$Pa$]', 'FontSize', smallFontSize, 'Interpreter', 'Latex');
                plot(xStorageMod(round(length(xStorageMod)/2):end),yStorageMod(round(length(xStorageMod)/2):end),'Linewidth',linewidthSize_loop)
                if CI
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1a = plot(xStorageMod(round(length(xStorageMod)/2):end),CI_storage(1,(round(length(xStorageMod)/2):end)),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1b = plot(xStorageMod(round(length(xStorageMod)/2):end),CI_storage(2,(round(length(xStorageMod)/2):end)),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                end
                set(gca,'xlim',[xStorageMod(round(length(xStorageMod)/2)) xStorageMod(end)])
                set(gca,'xticklabel',[])
    %             ylim auto
                hold off
                if i > 1
                    limSet = findobj(gca, '-property', 'ydata');
                    limSet = get(limSet, 'YData');
                    if iscell(limSet)
                        limSet = [limSet{:}];
                    end
                    set(gca,'ylim',[min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)])
                end

                subtightplot (2, 2, 3, [0.05 0.05], [0.125 0.05], [0.125 0.05]);
                hold on
                grid on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                set(findall(gcf,'-property','FontSize'),'FontSize',smallFontSize)
                set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop/2)
                xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', smallFontSize)
                ylabel('Loss Modulus [$Pa$]', 'FontSize', smallFontSize, 'Interpreter', 'Latex');
                plot(xLossMod(1:round(length(xLossMod)/2)),yLossMod(1:round(length(xLossMod)/2)),'Linewidth',linewidthSize_loop)
                if CI
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1a = plot(xLossMod(1:round(length(xLossMod)/2)),CI_loss(1,(1:round(length(xLossMod)/2))),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1b = plot(xLossMod(1:round(length(xLossMod)/2)),CI_loss(2,(1:round(length(xLossMod)/2))),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                end
                set(gca,'xlim',[xLossMod(1) xLossMod(round(length(xLossMod)/2))])
    %             ylim auto
                hold off
                if i > 1
                    limSet = findobj(gca, '-property', 'ydata');
                    limSet = get(limSet, 'YData');
                    if iscell(limSet)
                        limSet = [limSet{:}];
                    end
                    set(gca,'ylim',[min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)])
                end

                subtightplot (2, 2, 4, [0.05 0.05], [0.125 0.05], [0.15 0.025]);
                hold on
                grid on
                set(gca,'xscale','log')
                set(gca,'yscale','log')
                set(findall(gcf,'-property','FontSize'),'FontSize',smallFontSize)
                set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop/2)
                xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', smallFontSize)
    %             ylabel('Loss Modulus [$Pa$]', 'FontSize', smallFontSize, 'Interpreter', 'Latex');
                plot(xLossMod(round(length(xLossMod)/2):end),yLossMod(round(length(xLossMod)/2):end),'Linewidth',linewidthSize_loop)
                if CI
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1a = plot(xLossMod(round(length(xLossMod)/2):end),CI_loss(1,(round(length(xLossMod)/2):end)),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1a.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                    set(gca,'ColorOrderIndex',i) % Plot with the color we just used
                    h1b = plot(xLossMod(round(length(xLossMod)/2):end),CI_loss(2,(round(length(xLossMod)/2):end)),'Linewidth',linewidthSize_loop,'HandleVisibility','off');
                    h1b.Color(4) = fillAlpha; % UNDOCUMENTED: Change the opacity of the line we just plotted
                end
                set(gca,'xlim',[xLossMod(round(length(xLossMod)/2)) xLossMod(end)])
    %             ylim auto
                hold off
                if i > 1
                    limSet = findobj(gca, '-property', 'ydata');
                    limSet = get(limSet, 'YData');
                    if iscell(limSet)
                        limSet = [limSet{:}];
                    end
                    set(gca,'ylim',[min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)])
                end
            end
            
        case 'separate'
            figure(harmonicFig_storage)
            hold on
            if CI
                if isrow(xStorageMod)
                    fill([xStorageMod fliplr(xStorageMod)], [CI_storage(1,:) fliplr(CI_storage(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha,'HandleVisibility','off')
                else
                    fill([xStorageMod' fliplr(xStorageMod')], [CI_storage(1,:) fliplr(CI_storage(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha,'HandleVisibility','off')
                end
            end
            plot(xStorageMod,yStorageMod,'Linewidth',linewidthSize_loop,'Color',colorSpecList{i_dir})
            grid on
            set(gca,'xscale',xScaling)
            set(gca,'yscale',yScaling)
            set(gca,'Box',boxSetting)
            if i > 1
                limSet = findobj(gca, '-property', 'ydata');
                limSet = get(limSet, 'YData');
                if iscell(limSet)
                    limSet = [limSet{:}];
                end
                if CI
                   limSet = horzcat(limSet,real(CI_storage(1,:)),real(CI_storage(2,:)));
                end
                newLim = [min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)];
                currentLim = get(findobj(gca, '-property', 'ylim'),'YLim');
                if currentLim(1) > newLim(1) && newLim(1) > plotMin
                    set(gca,'ylim',[newLim(1) currentLim(2)]);
                    currentLim(1) = newLim(1);
                elseif currentLim(1) > newLim(1) && newLim(1) < plotMin
                    set(gca,'ylim',[plotMin currentLim(2)]);
                    currentLim(1) = plotMin;
                end
                if currentLim(2) < newLim(2)
                    set(gca,'ylim',[currentLim(1) newLim(2)]);
                    currentLim(2) = newLim(2);
                end
            end
            set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop)
            xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
            ylabel('Storage Modulus [$Pa$]', 'FontSize', bigFontSize, 'Interpreter', 'Latex');
            hold off

            figure(harmonicFig_loss)
            hold on
            if CI
                if isrow(xLossMod)
                    fill([xLossMod fliplr(xLossMod)], [CI_loss(1,:) fliplr(CI_loss(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha,'HandleVisibility','off')
                else
                    fill([xLossMod' fliplr(xLossMod')], [CI_loss(1,:) fliplr(CI_loss(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha,'HandleVisibility','off')
                end
            end
            plot(xLossMod,yLossMod,'Linewidth',linewidthSize_loop,'Color',colorSpecList{i_dir})
            grid on
            set(gca,'xscale',xScaling)
            set(gca,'yscale',yScaling)
            set(gca,'Box',boxSetting)
            if i > 1
                limSet = findobj(gca, '-property', 'ydata');
                limSet = get(limSet, 'YData');
                if iscell(limSet)
                    limSet = [limSet{:}];
                end
                if CI
                   limSet = horzcat(limSet,abs(CI_loss(1,:)),abs(CI_loss(2,:)));
                end
                newLim = [min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)];
                currentLim = get(findobj(gca, '-property', 'ylim'),'YLim');
                if currentLim(1) > newLim(1) && newLim(1) > plotMin
                    set(gca,'ylim',[newLim(1) currentLim(2)]);
                    currentLim(1) = newLim(1);
                elseif currentLim(1) > newLim(1) && newLim(1) < plotMin
                    set(gca,'ylim',[plotMin currentLim(2)]);
                    currentLim(1) = plotMin;
                end
                if currentLim(2) < newLim(2) && newLim(2) > plotMin
                    set(gca,'ylim',[currentLim(1) newLim(2)]);
                    currentLim(2) = newLim(2);
                end
            end
            set(findall(gcf,'-property','FontSize'),'FontSize',mediumFontSize)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',linewidthSize_loop)
            xlabel('Frequency [$Hz$]', 'Interpreter', 'latex', 'FontSize', bigFontSize)
            ylabel('Loss Modulus [$Pa$]', 'FontSize', bigFontSize, 'Interpreter', 'Latex');
            hold off
            
    end
    
    % Loss Angle Plot
    figure(lossAngleFig)
    if yyLossAngle
        if ismember(i_dir,leftPlotNums)
            yyaxis left
        elseif ismember(i_dir,rightPlotNums)
            yyaxis right
        end
    end
    hold on
    plot(xLossMod,yLossAngle,'Linewidth',linewidthSize,'linestyle','-','Color',colorSpecList{i_dir},'Marker','none')
    if exist('CI_angle','var')
        if isrow(xLossMod)
            fill([xLossMod fliplr(xLossMod)], [CI_angle(1,:) fliplr(CI_angle(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha, 'HandleVisibility', 'off', 'linestyle', '-', 'Marker', 'none')
        else
            fill([xLossMod' fliplr(xLossMod')], [CI_angle(1,:) fliplr(CI_angle(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha, 'HandleVisibility', 'off', 'linestyle', '-', 'Marker', 'none')
        end
    end
    set(gca,'Box',boxSetting)
    hold off
    
    % Force Plot
    % Figure 3, Subfigure b
    % Now, we create the force plot for our selected dataset      
    
    Files = dir([Folders{cellSelected(i_dir)} '\*dataStruct*.mat']);
    
    % Grab the info we will need below from the optimal set
    % Load the file
    clearvars dataStruct
    fileNo = optimalTermNumbers{cellSelected(i_dir)}(loadLevelSelected(i_dir));
    load([Files(fileNo).folder '\' Files(fileNo).name]);
    
    % Find the number of terms
    temp = strsplit(Files(fileNo).name,{'-','.'});
    termind = 3;
    if strcmpi(temp{1},'NIH')
        termind = termind+1;
    end
    temp = strsplit(temp{termind},{'_'});
    n_terms = str2num(temp{1});

    % Calculate the offsets
    indShift = 0;
    for k = 1:size(dataStruct,2)
        if isempty(dataStruct(k).t_average)
            indShift = indShift + 1;
        end
    end
    if strcmp(avAnalysis,'i')
        indShift = 0;
    end
    
    switch lower(modelSelected)
        case 'gm'
            if settingsStruct.smoothData == 'n'
                inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r,...
                    dataStruct(indShift+loadLevelSelected(i_dir)).h_r.^1.5];
                xDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).t_r;
                xDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).t_r;
                yDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).F_r;
                yDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).h_r;
            else
                inputs = [dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth,...
                    dataStruct(indShift+loadLevelSelected(i_dir)).h_r_smooth.^1.5];
                xDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth;
                xDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth;
                yDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).F_r_smooth;
                yDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).h_r_smooth;
            end
            
            timeVec = inputs(:,1);
            CI_force = NaN(2,length(timeVec));
            dt = dataStruct(indShift+loadLevelSelected(i_dir)).dt;
            st = timeVec(end);
            [~,~,h_conv,~,~] = ...
                makeGeneralizedMaxwellModel(fluidSetting_maxwell,elasticSetting_maxwell,...
                n_terms,timeVec,dt,st);
            modelParams = dataStruct(indShift+loadLevelSelected(i_dir)).maxwellParams;
            xModelForce = timeVec;
            yModelForce = alpha.*h_conv(modelParams,inputs);

            allParams = dataStruct(indShift+loadLevelSelected(i_dir)).maxwellParams_alternatives;
            allForce = NaN(length(timeVec),size(allParams,2));
            for j = 1:size(allParams,2)
                allForce(:,j) = alpha.*h_conv(allParams(:,j),inputs);
            end

            % Find the harmonic bounds represented by our datasets.
            CI_force(1,:) = max(allForce,[],2)';
            CI_force(2,:) = min(allForce,[],2)';
            CI_force = sort(CI_force,1,'descend');

        case 'gkv'
%                 allParams = dataStruct(indShift+loadLevelSelected(i_dir)).convParams_alternatives;
            % Do nothing. It's too labor intensive to get the maxwell
            % collocation results for every alternative.
            if settingsStruct.smoothData == 'n'
                inputs = [dataStruct(indShift+k).t_r,...
                    dataStruct(indShift+k).F_r];
                xDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).t_r;
                xDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).t_r;
                yDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).F_r;
                yDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).h_r;
            else
                inputs = [dataStruct(indShift+k).t_r_smooth,...
                    dataStruct(indShift+k).F_r_smooth];
                xDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth;
                xDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).t_r_smooth;
                yDataForce = dataStruct(indShift+loadLevelSelected(i_dir)).F_r_smooth;
                yDataInd = dataStruct(indShift+loadLevelSelected(i_dir)).h_r_smooth;
            end
            timeVec = inputs(:,1);
            CI_force = NaN(2,length(timeVec));
%                 CI_force(1,:) = NaN(size(timeVec));
%                 CI_force(2,:) = NaN(size(timeVec));
            xModelForce = timeVec;
            yModelForce = NaN(size(timeVec));

    end
    
    stepsize = floor(length(xModelForce)/10);
    
    figure(forceFig)
    hold on
    scatter(xDataForce,yDataForce,150,'*','LineWidth',linewidthSize/2,'MarkerEdgeColor',colorSpecList{i_dir})
    plot(xModelForce(1:stepsize:end),yModelForce(1:stepsize:end),'Linewidth',linewidthSize,'Color','k','Marker',markerSpecList{i_dir+1},'MarkerFaceColor',colorSpecList{i_dir},'MarkerSize',markerSize,'MarkerEdgeColor','black')
    fprintf('Force Data Figure, %s, %d Term(s), Standard Error: %4.4g\n',dataLabel{i_dir},n_terms,SError(yDataForce,yModelForce))

%     if plotHertz
%         plot(xModelHertz(1:stepsize:end),yModelHertz(1:stepsize:end),'--','Linewidth',linewidthSize,'Color',colorSpecList{i_dir})
%     end
    set(gca,'xscale','linear')
    set(gca,'yscale','linear')
    if CI
        if isrow(xModelForce)
            fill([xModelForce fliplr(xModelForce)], [CI_force(1,:) fliplr(CI_force(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha,'HandleVisibility','off')
        else
            fill([xModelForce' fliplr(xModelForce')], [CI_force(1,:) fliplr(CI_force(2,:))], colorSpecList{i_dir}, 'FaceAlpha', fillAlpha,'HandleVisibility','off')
        end
    end
    hold off
%     fprintf('Used %d Terms for the %s Force Plot.\n',n_terms,dataLabel{i_dir});
    
    figure(indFig)
    hold on
    scatter(xDataInd,yDataInd,150,'*','LineWidth',linewidthSize/2,'MarkerEdgeColor',colorSpecList{i_dir})
    set(gca,'xscale','linear')
    set(gca,'yscale','linear')
    hold off

    stepsize = floor(length(xModelHertz)/10);
    
    figure(hertzForceFig)
    hold on
    scatter(xDataHertz,yDataHertz,150,'*','LineWidth',linewidthSize/2,'MarkerEdgeColor',colorSpecList{i_dir})
    plot(xModelHertz(1:stepsize:end),yModelHertz(1:stepsize:end),'-','Linewidth',linewidthSize,'Color','k','Marker',markerSpecList{i_dir+1},'MarkerFaceColor',colorSpecList{i_dir},'MarkerSize',markerSize,'MarkerEdgeColor','black')
    fprintf('Hertz Force Figure, %s, Standard Error: %4.4g\n',dataLabel{i_dir},SError(yDataHertz,yModelHertz))
    set(gca,'xscale','linear')
    set(gca,'yscale','linear')
    hold off
    
end

externalLegend = {};
if ~isempty(plotExternal)

    externalData = cell(length(plotExternal));
    externalStiffness = cell(length(plotExternal));
    externalForce = cell(length(plotExternal));
    externalHarmonicLoss = cell(length(plotExternal));
    externalHarmonicStorage = cell(length(plotExternal));
    externalAbsMod = cell(length(plotExternal));
    externalLossAngle = cell(length(plotExternal));

    for i = 1:length(plotExternal)

        switch plotExternal{i}
            case 'Alcaraz2003'
                omega_alcaraz = 2*pi*logspace(-1,2,1000)';

                % A549
                G0 = 458; % Pa
                alpha = 0.22;
                mu = 1.68; % Pa-s
                omega_0 = 1; % s^-1
                eta = tan(alpha*pi/2);
                complexShear = G0.*(1+1i*eta).*(omega_alcaraz./omega_0).^alpha+1i.*omega_alcaraz.*mu;
                externalHarmonicLoss{i} = horzcat(externalHarmonicLoss{i},imag(complexShear));
                externalHarmonicStorage{i} = horzcat(externalHarmonicStorage{i},real(complexShear));
                externalLossAngle{i} = horzcat(externalLossAngle{i},atand(externalHarmonicLoss{i}./externalHarmonicStorage{i}));

                % Loss Angle Plot
                clearvars ax1_lossAng ax2_lossAng
                figure(lossAngleFig)
                axList = findobj(get(gcf, 'Children'), 'type', 'axes');
                ax1_lossAng = axList(1);
                numExpectedAxes = 1;
                if length(axList) == numExpectedAxes
                    if ~useLegend
                        ax2_lossAng = axes('Position',[0.65 0.65 0.225 0.25]);
                    else
                        ax2_lossAng = axes('Position',[0.25 0.65 0.125 0.25]);
                    end
                else
                    ax2_lossAng = axList(2);
                end
                currentColor = length(pathname)+1;
                box(ax2_lossAng, boxSetting)
                hold(ax2_lossAng, 'on')
                grid(ax2_lossAng, 'on')
                plot(ax2_lossAng,omega_alcaraz/(2*pi),externalLossAngle{i}(:,1),'Linewidth',linewidthSize*2, 'linestyle', '-', 'Color', colorOrder(currentColor,:))
                set(ax2_lossAng, 'xlim', [min(omega_alcaraz/(2*pi)) max(omega_alcaraz/(2*pi))])
                hold(ax2_lossAng, 'off')

                % BEAS-2B
                G0 = 496; % Pa
                alpha = 0.20;
                mu = 2.69; % Pa-s
                omega_0 = 1; % s^-1
                eta = tan(alpha*pi/2);
                complexShear = G0.*(1+1i*eta).*(omega_alcaraz./omega_0).^alpha+1i.*omega_alcaraz.*mu;
                externalHarmonicLoss{i} = horzcat(externalHarmonicLoss{i},imag(complexShear));
                externalHarmonicStorage{i} = horzcat(externalHarmonicStorage{i},real(complexShear));
                externalLossAngle{i} = horzcat(externalLossAngle{i},atand(externalHarmonicLoss{i}./externalHarmonicStorage{i}));

                % Loss Angle Plot
                hold(ax2_lossAng, 'on')
                plot(ax2_lossAng,omega_alcaraz/(2*pi),externalLossAngle{i}(:,2),'Linewidth',linewidthSize, 'linestyle', '-', 'Color', colorOrder(currentColor+1,:))
                hold(ax2_lossAng, 'off')

                set(findall(ax2_lossAng,'-property','FontSize'),'FontSize',smallFontSize,...
                    'FontWeight','bold')
                set(ax2_lossAng,'xscale',xScaling)
                set(ax2_lossAng,'TickLength',[0.02 0.02]/2,'LineWidth',linewidthSize/2)
                
                externalLegend = horzcat(externalLegend, 'A549 (Alcaraz 2003)', 'BEAS-2B (Alcaraz 2003)');

        end

    end
end

% Create the Hertzian box and whisker plot
figure(hertzFig)
hold on
h = boxplot(hertzList./1e3,hertzLabel,'OutlierSize',markerSize,'Symbol','rx');
set(h,{'linewidth'},{4})
h = gca;
h.XAxis.TickLabelInterpreter = 'latex';
h = findobj(hertzFig,'tag','Outliers');
set(gca,'xscale','linear')
set(gca,'Box',boxSetting)
set(gca,'yscale','linear')
hertzLabelNum = ones(size(hertzLabel));
hertzLabelNum(strcmpi(hertzLabel,'Melanomas')) = 2;
hertzLabelNum(strcmpi(hertzLabel,'Fibroblasts')) = 3;

outData = flip(get(h,'YData'));
valsToRemove = ismember(hertzList./1e3,(horzcat(outData{1},outData{2},outData{3})));
scatter(hertzLabelNum(~valsToRemove & hertzLabelNum==1),hertzList(~valsToRemove & hertzLabelNum==1)./1e3,markerSize*10,[colorSpecList{1} 'o'],'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',1,'MarkerEdgeColor',colorSpecList{1},'jitter','on','jitterAmount',0.1)
scatter(hertzLabelNum(~valsToRemove & hertzLabelNum==2),hertzList(~valsToRemove & hertzLabelNum==2)./1e3,markerSize*10,[colorSpecList{2} 'o'],'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',1,'MarkerEdgeColor',colorSpecList{2},'jitter','on','jitterAmount',0.1)
scatter(hertzLabelNum(~valsToRemove & hertzLabelNum==3),hertzList(~valsToRemove & hertzLabelNum==3)./1e3,markerSize*10,[colorSpecList{3} 'o'],'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',1,'MarkerEdgeColor',colorSpecList{3},'jitter','on','jitterAmount',0.1)

cellMean = mean(hertzList(strcmpi(hertzLabel,'Melanocytes'))./1e3);
cellSD = std(hertzList(strcmpi(hertzLabel,'Melanocytes'))./1e3);
% text(1,-0.55,sprintf('%1.3f\\pm%1.3fkPa',cellMean,cellSD),'HorizontalAlignment','center','FontWeight','normal','FontSize',smallFontSize);

cellMean = mean(hertzList(strcmpi(hertzLabel,'Melanomas'))./1e3);
cellSD = std(hertzList(strcmpi(hertzLabel,'Melanomas'))./1e3);
% text(2,-0.55,sprintf('%1.3f\\pm%1.3f kPa',cellMean,cellSD),'HorizontalAlignment','center','FontWeight','normal','FontSize',smallFontSize);

cellMean = mean(hertzList(strcmpi(hertzLabel,'Fibroblasts'))./1e3);
cellSD = std(hertzList(strcmpi(hertzLabel,'Fibroblasts'))./1e3);
% text(3,-0.55,sprintf('%1.3f\\pm%1.3f kPa',cellMean,cellSD),'HorizontalAlignment','center','FontWeight','normal','FontSize',smallFontSize);

yt = get(gca, 'YTick');
axis([0.5 3.5 0 ceil(max(hertzList./1e3))+1.3])
lineOffset = 0.00925;

[~,cellP] = ttest2(hertzList(strcmpi(hertzLabel,'Melanocytes'))./1e3,hertzList(strcmpi(hertzLabel,'Melanomas'))./1e3,'Tail','both','Alpha',0.05,'Vartype','unequal');
if cellP < 0.0001
    pString = sprintf('p < 0.0001');
else
    pString = sprintf('p = %1.4f',cellP);
end
plot([1 2],[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Melanomas'}))./1e3)),'b-','Linewidth',4)
plot([1+lineOffset 1+lineOffset],[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Melanomas'}))./1e3))-[0 0.1],'b-','Linewidth',4)
plot([2-lineOffset 2-lineOffset],[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Melanomas'}))./1e3))-[0 0.1],'b-','Linewidth',4)
text(1.5,0.25+ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Melanomas'}))./1e3)),pString,'HorizontalAlignment','center','FontWeight','bold','FontSize',smallFontSize);

[~,cellP] = ttest2(hertzList(strcmpi(hertzLabel,'Melanocytes'))./1e3,hertzList(strcmpi(hertzLabel,'Fibroblasts'))./1e3,'Tail','both','Alpha',0.05,'Vartype','unequal');
if cellP < 0.0001
    pString = sprintf('p < 0.0001');
else
    pString = sprintf('p = %1.4f',cellP);
end
plot([1 3],0.6+[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Fibroblasts'}))./1e3)),'b-','Linewidth',4)
plot([1+lineOffset 1+lineOffset],0.6+[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Fibroblasts'}))./1e3))-[0 0.1],'b-','Linewidth',4)
plot([3-lineOffset 3-lineOffset],0.6+[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Fibroblasts'}))./1e3))-[0 0.1],'b-','Linewidth',4)
text(2,0.85+ceil(max(hertzList(ismember(hertzLabel,{'Melanocytes','Fibroblasts'}))./1e3)),pString,'HorizontalAlignment','center','FontWeight','bold','FontSize',smallFontSize);

[~,cellP] = ttest2(hertzList(strcmpi(hertzLabel,'Melanomas'))./1e3,hertzList(strcmpi(hertzLabel,'Fibroblasts'))./1e3,'Tail','both','Alpha',0.05,'Vartype','unequal');
if cellP < 0.0001
    pString = sprintf('p < 0.0001');
else
    pString = sprintf('p = %1.4f',cellP);
end
plot([2 3],[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanomas','Fibroblasts'}))./1e3)),'b-','Linewidth',4)
plot([2+lineOffset 2+lineOffset],[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanomas','Fibroblasts'}))./1e3))-[0 0.1],'b-','Linewidth',4)
plot([3-lineOffset 3-lineOffset],[1 1]*ceil(max(hertzList(ismember(hertzLabel,{'Melanomas','Fibroblasts'}))./1e3))-[0 0.1],'b-','Linewidth',4)
text(2.5,0.25+ceil(max(hertzList(ismember(hertzLabel,{'Melanomas','Fibroblasts'}))./1e3)),pString,'HorizontalAlignment','center','FontWeight','bold','FontSize',smallFontSize);

hold off

%% Save the plots
% Save Harmonic Plot
if ~strcmp(harmonicConfig,'separate')
    figure(harmonicFig)
    if useLegend
%         if stackHarmonics
%             legend(horzcat(legendString{2:2:end},externalLegend),'location','eastoutside','orientation','horizontal','numColumns',1)
%         else
            legend(horzcat(horzcat(legendString{2:2:end},externalLegend),horzcat(legendString{2:2:end},externalLegend)),'location','eastoutside','orientation','horizontal','numColumns',1)
%         end    
    end
    saveas(harmonicFig, [outputpathname '\' sprintf('%s-Harmonics',saveLabel)], 'jpg')
    saveas(harmonicFig, [outputpathname '\' sprintf('%s-Harmonics',saveLabel)], 'fig')
else
    figure(harmonicFig_storage)
%     limSet = findobj(gca, '-property', 'ydata');
%     limSet = get(limSet, 'YData');
%     limSet = [limSet{:}];
%     set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
    limSet = findobj(gca, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
    if useLegend
%         if stackHarmonics
%             legend(horzcat(legendString{2:2:end},externalLegend),'location','eastoutside','orientation','horizontal','numColumns',1)
%         else
            legend(horzcat(horzcat(legendString{2:2:end},externalLegend),horzcat(legendString{2:2:end},externalLegend)),'location','eastoutside','orientation','horizontal','numColumns',1)
%         end    
    end
    saveas(harmonicFig_storage, [outputpathname '\' sprintf('%s-StorageModulus',saveLabel)], 'jpg')
    saveas(harmonicFig_storage, [outputpathname '\' sprintf('%s-StorageModulus',saveLabel)], 'fig')
    
    figure(harmonicFig_loss)
%     limSet = findobj(gca, '-property', 'ydata');
%     limSet = get(limSet, 'YData');
%     limSet = [limSet{:}];
% 	set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
    limSet = findobj(gca, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
    if useLegend
%         if stackHarmonics
%             legend(horzcat(legendString{2:2:end},externalLegend),'location','eastoutside','orientation','horizontal','numColumns',1)
%         else
            legend(horzcat(horzcat(legendString{2:2:end},externalLegend),horzcat(legendString{2:2:end},externalLegend)),'location','eastoutside','orientation','horizontal','numColumns',1)
%         end    
    end
    saveas(harmonicFig_loss, [outputpathname '\' sprintf('%s-LossModulus',saveLabel)], 'jpg')
    saveas(harmonicFig_loss, [outputpathname '\' sprintf('%s-LossModulus',saveLabel)], 'fig')
end

% Save Loss Angle Plot
figure(lossAngleFig)
axList = (findobj(get(gcf, 'Children'), '-depth', 1, 'type', 'axes'));
if length(axList) > 1
    mainAx = axList(2);
    subAx = axList(1);
    if yyLossAngle
        yyaxis(mainAx, 'left')
        temp1 = flipud(mainAx.Children);
        yyaxis(mainAx, 'right')
        temp2 = flipud(mainAx.Children);
        mainAxLeg = vertcat(temp1,temp2);
        mainAxOrder = horzcat(leftPlotNums,rightPlotNums);
        mainAxLeg = mainAxLeg(mainAxOrder);
        
        yyaxis(mainAx, 'left')
        limSet = findobj(gca, '-property', 'ydata');
        limSet = get(limSet, 'YData');
        limSet = [limSet{:}];
        set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
        yyaxis(mainAx, 'right')
        limSet = findobj(gca, '-property', 'ydata');
        limSet = get(limSet, 'YData');
        limSet = [limSet{:}];
        set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
        
    else
        mainAxLeg = flipud(findall(get(mainAx, 'Children'), 'type', 'line'));
        
        limSet = findobj(mainAx, '-property', 'ydata');
        limSet = get(limSet, 'YData');
        limSet = [limSet{:}];
        set(mainAx,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
    end
    legendObjs = vertcat(mainAxLeg,...
        flipud(findall(get(subAx, 'Children'), 'type', 'line')));
    
    limSet = findobj(mainAx, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(mainAx,'xtick',sort(10.^unique(floor(log10(limSet)))));
    
    limSet = findobj(subAx, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(subAx,'xtick',sort(10.^unique(floor(log10(limSet)))));
else
    mainAx = axList;
    legendObjs = flipud(findall(get(mainAx, 'Children'), 'type', 'line'));
    limSet = findobj(mainAx, '-property', 'ydata');
    limSet = get(limSet, 'YData');
    limSet = [limSet{:}];
    set(mainAx,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
    limSet = findobj(mainAx, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(mainAx,'xtick',sort(10.^unique(floor(log10(limSet)))));
end
if useLegend
    legend(legendObjs,horzcat(legendString,externalLegend),'location','eastoutside','orientation','vertical','numColumns',1)
end
if ~useLegend
    legendFig = figure;
    legendEntries = horzcat(legendString,externalLegend);
    hold on
    % Plot whatever you like
    x = 1:10;
    y = NaN;
    for i = 1:length(Files)
        plot(x, y .* x, 'DisplayName', legendEntries{i},'color',legendObjs(i).Color,'lineWidth',legendObjs(i).LineWidth,...
            'LineStyle',legendObjs(i).LineStyle,'Marker',legendObjs(i).Marker,'MarkerSize',legendObjs(i).MarkerSize,...
            'MarkerFaceColor',legendObjs(i).MarkerFaceColor)
    end
    % Initial values to capture the entire legend
    % Should fit most modern screens
    set(gcf,'Position',[0,0,1024,1024]);
    % Call the legend to your choice, I used a horizontal legend here
    legend_handle = legend('Orientation','vertical','FontSize',mediumFontSize);
    % Set the figure Position using the normalized legend Position vector
    % as a multiplier to the figure's current position in pixels
    % This sets the figure to have the same size as the legend
    set(gcf,'Position',(get(legend_handle,'Position')...
        .*[0, 0, 1, 1].*get(gcf,'Position')));
    % The legend is still offset so set its normalized position vector to
    % fill the figure
    set(legend_handle,'Position',[0,0,1,1]);
    % Put the figure back in the middle screen area
    set(gcf, 'Position', get(gcf,'Position') + [500, 400, 0, 0]);

    saveas(legendFig, [outputpathname '\' sprintf('%s-LossAngleLegend',saveLabel)], 'jpg')
    saveas(legendFig, [outputpathname '\' sprintf('%s-LossAngleLegend',saveLabel)], 'fig')
    close(legendFig)
    clearvars legendFig
end
saveas(lossAngleFig, [outputpathname '\' sprintf('%s-LossAngle',saveLabel)], 'jpg')
saveas(lossAngleFig, [outputpathname '\' sprintf('%s-LossAngle',saveLabel)], 'fig')

% Save Force Plot
figure(forceFig)
limSet = findobj(gca, '-property', 'ydata');
limSet = get(limSet, 'YData');
limSet = [limSet{:}];
set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
if strcmp(get(findobj(gca, '-property', 'xscale'),'XScale'),'log')
    limSet = findobj(gca, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
end
if useLegend
    legend(legendString,'location','eastoutside','orientation','horizontal','numColumns',1)
end
saveas(forceFig, [outputpathname '\' sprintf('%s-ForcePlot',saveLabel)], 'jpg')
saveas(forceFig, [outputpathname '\' sprintf('%s-ForcePlot',saveLabel)], 'fig')

% Save Indentation Plot
figure(indFig)
limSet = findobj(gca, '-property', 'ydata');
limSet = get(limSet, 'YData');
limSet = [limSet{:}];
set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
if strcmp(get(findobj(gca, '-property', 'xscale'),'XScale'),'log')
    limSet = findobj(gca, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
end
if useLegend
    legend(legendString,'location','eastoutside','orientation','horizontal','numColumns',1)
end
saveas(indFig, [outputpathname '\' sprintf('%s-IndPlot',saveLabel)], 'jpg')
saveas(indFig, [outputpathname '\' sprintf('%s-IndPlot',saveLabel)], 'fig')

% Save Hertz Plot
figure(hertzFig)
saveas(hertzFig, [outputpathname '\' sprintf('%s-HertzPlot',saveLabel)], 'jpg')
saveas(hertzFig, [outputpathname '\' sprintf('%s-HertzPlot',saveLabel)], 'fig')

% Save Hertz Force Plot
figure(hertzForceFig)
limSet = findobj(gca, '-property', 'ydata');
limSet = get(limSet, 'YData');
limSet = [limSet{:}];
set(gca,'ylim',sort([min(limSet,[],'all')*(1-limMargin) max(limSet,[],'all')*(1+limMargin)]))
if strcmp(get(findobj(gca, '-property', 'xscale'),'XScale'),'log')
    limSet = findobj(gca, '-property', 'xdata');
    limSet = get(limSet, 'XData');
    limSet = [limSet{:}];
    set(gca,'xtick',sort(10.^unique(floor(log10(limSet)))));
end
if useLegend
    legend(legendString,'location','eastoutside','orientation','horizontal','numColumns',1)
end
saveas(hertzForceFig, [outputpathname '\' sprintf('%s-HertzForcePlot',saveLabel)], 'jpg')
saveas(hertzForceFig, [outputpathname '\' sprintf('%s-HertzForcePlot',saveLabel)], 'fig')

close all

%% Open Initially Selected Directory
clearvars -except outputpathname
close all
disp('Finished Plotting the Paper Figures.')
winopen(outputpathname);