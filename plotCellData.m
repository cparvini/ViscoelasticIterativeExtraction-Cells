clear all
close all
clc

% ======================== %
% Cell Data Visualization  %
%          Script          %  
% ======================== %
% This script is used to 
% visualize the cell data
% contained in the data 
% directory.
% ======================== %

addpath(genpath(pwd))

% Pick the AFM Data Directory and Choose Data Extraction Settings
originalPath = uigetdir(pwd,...
        'Select the Folder Containing Your AFM Files');
    
% outputpathname = uigetdir(originalPath,...
%         'Select an Output Folder');

% Check to see if there are subdirectories
dirContents = dir(originalPath);
subFolders = dirContents([dirContents.isdir]);
subFolders(contains({subFolders.name}, {'.','..','Plots'})) = [];

% If the user provides a main directory with many subdirectories containing
% data, we should loop through all directories and analyze each in turn.
if ~isempty(subFolders)
    Folders = cell(1,length(subFolders));
    Folders = cellfun(@(root,sub)[root '\' sub],{subFolders.folder},{subFolders.name},'UniformOutput',false);
else
    Folders = {originalPath};
end

for i_dir = 1:length(Folders)
    
    % Start with a clean slate
    clearvars -except i_dir Folders originalPath
    close all
    clc
    fprintf('Analyzing Directory #%d of %d\n',i_dir,length(Folders));
    
    % Use the current 
    path = Folders{i_dir};
    
    % Alternatively, set those values manually
    minTimescale = 1e-4;                    % This is the time value on which the first viscoelastic element will be centered
    tipGeom = "spherical";                  % The experiment tip geometry for the files that are being loaded
    
    % Settings for how to use the loaded data during analysis
    useSmoothData = 0;                      % The user can choose to use filtered data, or the original raw data for fitting
    useAveragedData = 1;                    % Choose whether to average force curves with the same approach velocity 
    
    % Settings for loading the data
    loadDataSettings  = struct();
    
    % Required Settings:
    loadDataSettings.includeRetract = 1;             % Don't include data from the retract curve
    loadDataSettings.filterType = 'none';            % Choose the filter used to smooth data
    loadDataSettings.findRep = 'forward';            % Search direction for the repulsive region
    loadDataSettings.removeNegatives = 1;            % Remove negative values in the data stream
    
    % Conditional Settings (depending on filter):
    loadDataSettings.N = 2;                          % Order of Butterworth filter (if used)
    loadDataSettings.cutoff_Hz = 5000;               % Cutoff frequency of Butterworth (if used)
    
    % Load the AFM Data
    dataStruct = LoadAFMData(path,loadDataSettings);
    
    % Choose the right data from our structure
    % To begin, we need to either select the averaged data rows (stored at the
    % bottom of the structure), and know the number of files used to make them
    % (the rows leading up to the averages). There will be one averaged
    % datasets per approach velocity.
    numFiles = 0;
    avgCount = 0;
    for i = 1:size(dataStruct,2)
        if isempty(dataStruct(i).t_average)
            numFiles = numFiles + 1;
        else
            avgCount = avgCount + 1;
        end
        v_approach(i) = dataStruct(i).v_approach;
    end

    % Get our offsets, knowing how many files and averages there are
    if useAveragedData
        indShift = numFiles;
        loopMax = avgCount;
    else
        indShift = 0;
        loopMax = numFiles;
    end
    
    % Get our unique velocities
    v_approach = round(v_approach,2,'significant');
    v_unique = (unique(v_approach));
    
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
        legend(legendString,'Location','eastoutside','Orientation','Vertical', 'FontSize', 10)
        hold off

        saveas(gcf, [path sprintf('/AFM-SFS-CorrectedDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [path sprintf('/AFM-SFS-CorrectedDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end

    clearvars legendString titleString

    for i = 1:length(v_unique)

        legendString = {};
        figure(2)
        clf
        set(gcf, 'Position', get(0, 'Screensize'));
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
                else
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

        saveas(gcf, [path sprintf('/AFM-SFS-RepulsiveDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [path sprintf('/AFM-SFS-RepulsiveDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end

    for i = 1:length(v_unique)

        legendString = {};
        figure(3)
        clf
        set(gcf, 'Position', get(0, 'Screensize'));
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
                else
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

        saveas(gcf, [path sprintf('/AFM-SFS-IndentationDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [path sprintf('/AFM-SFS-IndentationDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end

    for i = 1:length(v_unique)

        legendString = {};
        figure(4)
        clf
        set(gcf, 'Position', get(0, 'Screensize'));
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
                else
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

        saveas(gcf, [path sprintf('/AFM-SFS-ForceDataset-Velocity_%.f_nm-s',v_unique(i)/1E-9)], 'jpg')
        savefig(gcf, [path sprintf('/AFM-SFS-ForceDataset-Velocity_%.f_nm-s.fig',v_unique(i)/1E-9)],'compact')
    end
    
end