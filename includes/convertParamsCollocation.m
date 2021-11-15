function [beta] = convertParamsCollocation(inputParams,outputModel,elasticSetting_input,fluidSetting_input,elasticSetting_output,fluidSetting_output,n_terms,padSize,timeVec,dt,st,minTimescale,forwardFitTimescale)
%convertParams Determine Conjugate Model Parameters by Fitting Harmonics
%   This function is used to convert input parameters from one model into
%   the corresponding parameters for another. It is designed for use only
%   with the Generalized Voigt and Generalized Maxwell models, but is also
%   useful for the Schapery-Framework models (since they are either based
%   on Voigt or Maxwell models). To assess the quality of fit, set
%   "showPlots" to 1 such that figures are created for the storage and loss
%   moduli.

lsqoptions = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective',...
    'MaxFunctionEvaluations', 1e5,...
    'MaxIterations', 1e5,...
    'FiniteDifferenceType','central',...
    'FunctionTolerance', 1e-60,...
    'OptimalityTolerance', 1e-60,...
    'StepTolerance', 1e-60,...
    'Display', 'none');

opts = optimoptions('fmincon','Algorithm','interior-point',...
    'MaxFunctionEvaluations', 5e4,...
    'MaxIterations', 5e4,...
    'Display','none',...
    'FiniteDifferenceType','central',...
    'FunctionTolerance', 1e-60,...
    'OptimalityTolerance',1e-60,...
    'StepTolerance', 1e-60);

showPlots = 1;
maxNum = 1e2;
linewidthSize = 5;
markerSize = linewidthSize*25;
loopType = 'parfor'; % Choose 'while' or 'parfor'

% Make array of log-spaced S values
% a = floor(log10(1/(10^(ceil(log10(timeVec(end)))))))-10;
% b = ceil(log10(1/(10^(floor(log10(dt))))))+5;
a = -100;
b = 100;
% s = sigma + i*omega
% sig = logspace(a,b,(b-a)*10);
sig = 0;
wi = 1i*( 2*pi*logspace(a,b,(b-a)*5) );
slog = sig + wi;
% slog = logspace(a,b);
% [sigmesh,wimesh] = meshgrid(sig,wi);
% slog = sigmesh+wimesh;

switch lower(outputModel)
    case 'voigt'
        % The input parameters are stiffnesses
        [~,~,~,lb,ub,~,~] = makeGeneralizedVoigtModel(elasticSetting_output,fluidSetting_output,n_terms,timeVec,dt,st,minTimescale,forwardFitTimescale);
        
        if strcmp(elasticSetting_output,'y')
            if strcmp(fluidSetting_output,'y')
                tauInds = (3:2:length(ub)-1);
                modulusInds = horzcat(1,(2:2:length(ub)-1));
            else
                tauInds = (3:2:length(ub));
                modulusInds = horzcat(1,(2:2:length(ub)));
            end
        else
            if strcmp(fluidSetting_output,'y')
                tauInds = (2:2:length(ub)-1);
                modulusInds = (1:2:length(ub)-1);
            else
                tauInds = (2:2:length(ub));
                modulusInds = (1:2:length(ub));
            end
        end
        
        % Free up the relaxation times for proper adjustment
        lb(tauInds) = 10^(floor(log10(dt))-1);
        ub(tauInds) = 10^(ceil(log10(timeVec(end)))+1);
        
        % The moduli should still be reasonable, so leave those indices
        % alone
        
        % Make U(s)
        UArms = cell(1,ceil(n_terms + padSize));
        for iii = 1:size(UArms,2)
            if iii == 1 && elasticSetting_output == 'y'
                UArms{iii} = sprintf('c(1)'); % E_e
            elseif iii == size(UArms,2) && fluidSetting_output == 'y'
                UArms{iii} = sprintf('0'); % phi
            else
                if strcmp(elasticSetting_output, 'y')
                    UArms{iii} = sprintf('(c(%d)) ./ (1.0 + ( c(%d).*s ) )',(iii-1)*2,(iii-1)*2+1);
                else
                    UArms{iii} = sprintf('(c(%d)) ./ (1.0 + ( c(%d).*s ) )',(iii)*2-1,(iii)*2);
                end
            end
        end

        % Create the Maxwell storage element series
        tempString1 = '(';
        for iii = 1:size(UArms,2)
            if iii < size(UArms,2)
                tempString1 = horzcat(tempString1, sprintf('(%s)+',UArms{iii}));
            else
                tempString1 = horzcat(tempString1, sprintf('(%s))',UArms{iii}));
            end
        end
        
        UString = sprintf('@(s,c) %s',tempString1);
        Ufunc = str2func(UString);
        
        % Make Q(s)
        QArms = cell(1,ceil(n_terms + padSize));
        for iii = 1:size(QArms,2)
            if iii == 1 && elasticSetting_input == 'y'
                QArms{iii} = sprintf('c(1)'); % E_e
            elseif iii == size(QArms,2) && fluidSetting_input == 'y'
                QArms{iii} = sprintf('0'); % phi
            else
                if strcmp(elasticSetting_input, 'y')
                    QArms{iii} = sprintf('(c(%d).*c(%d).*s) ./ (1.0 + ( c(%d).*s ) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
                else
                    QArms{iii} = sprintf('(c(%d).*c(%d).*s) ./ (1.0 + ( c(%d).*s ) )',(iii)*2-1,(iii)*2,(iii)*2);
                end
            end
        end

        % Create the Maxwell storage element series
        tempString1 = '(';
        for iii = 1:size(QArms,2)
            if iii < size(QArms,2)
                tempString1 = horzcat(tempString1, sprintf('(%s)+',QArms{iii}));
            else
                tempString1 = horzcat(tempString1, sprintf('(%s))',QArms{iii}));
            end
        end
        
        QString = sprintf('@(s,c) %s',tempString1);
        Qfunc = str2func(QString);
        
        % Create Y data and Models
        y_data = Qfunc(slog,inputParams);
        outFunc = @(s,c) 1./Ufunc(s,c);
        
        eflag = 0;
        
        beta0 = NaN(size(ub));
%         beta0(tauInds) = 10.^(log10(lb(tauInds))+rand(size(tauInds)).*(log10(ub(tauInds))-log10(lb(tauInds))));
%         beta0(modulusInds) = 10.^(log10(lb(modulusInds))+rand(size(modulusInds)).*(log10(ub(modulusInds))-log10(lb(modulusInds))));
%         if modulusInds(1) == 1
%             beta0(1) = 1./(inputParams(1));
%         end
        beta0(tauInds) = inputParams(tauInds);
        beta0(modulusInds) = 10.^(log10(1./(10.^log10(abs(max(y_data,[],'all')))))+rand(size(modulusInds)).*(log10(1./(10.^log10(abs(min(y_data,[],'all')))))-log10(1./(10.^log10(abs(max(y_data,[],'all')))))));
        
        if strcmp(loopType,'while')
            runNum = 0;
            bestRes = Inf;
            bestParams = NaN(size(ub));
        
            while any(eflag == [0 2 3])
                runNum = runNum + 1;

                objFunc = @(c,s) vertcat(real(outFunc(s,c)),imag(outFunc(s,c)));
                [beta,resnorm,~,eflag,~,~,~] = ...
                    lsqcurvefit(objFunc,beta0,slog,vertcat(real(y_data),imag(y_data)),lb,[],lsqoptions);

    %             costFunc = @(c) sum(abs((outFunc(slog,c)-y_data)./y_data),'all');
    %             [beta,~,eflag,~] = ...
    %                 fmincon(costFunc,beta0,[],[],[],[],lb,[],[],opts);
    %             resnorm = sum(abs((outFunc(slog,beta)-y_data)./y_data),'all');

                % Show the performance
    %             figure
    %             plot(log10(slog),log10(y_data),'*')
    %             hold on
    %             plot(log10(slog),log10(outFunc(slog,beta0)),'b')
    %             plot(log10(slog),log10(outFunc(slog,beta)),'r')
    %             hold off

                beta0 = NaN(size(ub));
    %             beta0(tauInds) = 10.^(log10(lb(tauInds))+rand(size(tauInds)).*(log10(ub(tauInds))-log10(lb(tauInds))));
    %             beta0(modulusInds) = 10.^(log10(lb(modulusInds))+rand(size(modulusInds)).*(log10(ub(modulusInds))-log10(lb(modulusInds))));
                beta0(tauInds) = inputParams(tauInds);
                beta0(modulusInds) = 10.^(log10(1./(10.^log10(abs(max(y_data,[],'all')))))+rand(size(modulusInds)).*(log10(1./(10.^log10(abs(min(y_data,[],'all')))))-log10(1./(10.^log10(abs(max(y_data,[],'all')))))));

                if resnorm < bestRes
                    bestRes = resnorm;
                    bestParams = beta;
                    if showPlots
                        fprintf('New Best Params!\n');
                    end
                end

                if showPlots
                    fprintf('Iteration %d; Exit Flag %d\n',runNum,eflag);
                end

                if runNum >= maxNum
                    if showPlots
                        fprintf('\nMaximum No. of Translation Attempts Reached...The New Predicted Parameters May Be Inaccurate.\n');
                    end
                    break;
                end
            end
        else
            
            if isempty(gcp('nocreate'))
                %parpool(N_workers)
                parpool(8,'IdleTimeout', Inf)
            end
            
            beta_dist = zeros(length(ub),maxNum);
            beta0_dist = zeros(size(beta_dist));
            resnorm_dist = zeros(1,maxNum);
            objFunc = @(c,s) vertcat(real(outFunc(s,c)),imag(outFunc(s,c)));
            
            progressString = sprintf('Converting GM to GKV\nInvestigating S-Plane Functions\nParallel Search Running...');
            hbar = parfor_progressbar(maxNum,progressString);
            warning('off');
            
            parfor i = 1:maxNum
                
                beta0 = NaN(size(ub));
                beta0(tauInds) = 10.^(log10(lb(tauInds))+rand(size(tauInds)).*(log10(ub(tauInds))-log10(lb(tauInds))));
                beta0(modulusInds) = 10.^(log10(lb(modulusInds))+rand(size(modulusInds)).*(log10(ub(modulusInds))-log10(lb(modulusInds))));
%                 beta0(tauInds) = inputParams(tauInds);
%                 beta0(modulusInds) = 10.^(log10(1./(10.^log10(abs(max(y_data,[],'all')))))+rand(size(modulusInds)).*(log10(1./(10.^log10(abs(min(y_data,[],'all')))))-log10(1./(10.^log10(abs(max(y_data,[],'all')))))));
%                 if modulusInds(1) == 1
%                     beta0(1) = 1./(inputParams(1));
%                 end
                
                [beta,resnorm,~,~,~,~,~] = ...
                    lsqcurvefit(objFunc,beta0,slog,vertcat(real(y_data),imag(y_data)),lb,[],lsqoptions);
                
                beta_dist(:,i) = beta;
                beta0_dist(:,i) = beta0;
                resnorm_dist(i) = resnorm;
                hbar.iterate(1) % Increase progressbar by 1 iteration
            end
            
            close(hbar)
            warning('on');
            
            [~,idx] = min(resnorm_dist(resnorm_dist > 0));
            bestBeta0 = beta0_dist(:,idx);
            bestParams = beta_dist(:,idx);
            
        end
        
        beta = bestParams;
        
        if showPlots
            % Show the performance
            fittingResultTest = figure('position',[500 20 1000 1000]);
            plot(log10(slog),log10(y_data),'*')
            hold on
            if strcmp(loopType,'parfor')
%                 plot(log10(slog),log10(outFunc(slog,bestBeta0)),'b')
            end
            plot(log10(slog),log10(outFunc(slog,beta)),'r')
            title('Fitting Result: Conversion to GM', 'FontSize', 16)
            set(findall(gcf,'-property','FontSize'),'FontSize',34)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',5)
            xlabel('$s, log(1/s)$', 'FontSize', 40, 'Interpreter', 'Latex')
            ylabel('$U(s), log(1/Pa)$', 'FontSize', 40, 'Interpreter', 'Latex')
            hold off
            
            % First, generate a frequency array in log scale
            de0 = 2.0.*pi.*(1./timeVec(end));
            maxi = 2.0.*pi.*(1./dt);
            omega = log_tw(de0,maxi);
            
            Jstorage = J_storage_advanced(omega,beta,elasticSetting_output,fluidSetting_output);
            Jloss = J_loss_advanced(omega,beta,elasticSetting_output,fluidSetting_output);
            Jabs = sqrt(Jstorage.^2+Jloss.^2);
            
            [~,~,~,~,padSize] = makeGeneralizedMaxwellModel(fluidSetting_input,elasticSetting_input,n_terms,timeVec,dt,st);
            [E_storage,E_loss,~,~,~,~] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_input,elasticSetting_input,n_terms,timeVec,minTimescale,padSize,forwardFitTimescale);
            
            if exist('storagePlotTest','var')
                clearvars storagePlotTest
            end
            if exist('lossModulusPlotTest','var')
                clearvars lossModulusPlot
            end
            storagePlotTest = figure('position',[500 30 1000 1000]);
            hold on
            grid on
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            title('Storage Modulus: Conversion to GKV', 'FontSize', 16)
            set(findall(gcf,'-property','FontSize'),'FontSize',34)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',5)
            xlabel('$\omega, \,rad/s$', 'FontSize', 40, 'Interpreter', 'Latex')
            ylabel('$E''(\omega)$', 'FontSize', 40, 'Interpreter', 'Latex')
            % xlim([1./10^ceil(log10(dataStruct(numFiles+kk).t_r(end))) 1./10^ceil(log10(dataStruct(numFiles+kk).dt))])
            hold off

            lossModulusPlotTest = figure('position',[500 40 1000 1000]);
            hold on
            grid on
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            title('Loss Modulus: Conversion to GKV', 'FontSize', 16)
            set(findall(gcf,'-property','FontSize'),'FontSize',34)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',5)
            xlabel('$\omega, \,rad/s$', 'FontSize', 40, 'Interpreter', 'Latex')
            ylabel('$E''''(\omega)$', 'FontSize', 40, 'Interpreter', 'Latex')
            % xlim([1./10^ceil(log10(dataStruct(numFiles+kk).t_r(end))) 1./10^ceil(log10(dataStruct(numFiles+kk).dt))])
            hold off
            
            % Calculate E_storage
            figure(storagePlotTest);
            hold on
            plot(omega, Jstorage./(Jabs.^2), '-r', 'LineWidth', 7)
            scatter(omega, E_storage(omega,beta),markerSize,'o','Linewidth',linewidthSize/2)
            hold off

            % Calculate E_loss
            figure(lossModulusPlotTest);
            hold on
            plot(omega, Jloss./(Jabs.^2), '-r', 'LineWidth', 7)
            scatter(omega, E_loss(omega,beta),markerSize,'o','Linewidth',linewidthSize/2)
            hold off
        end
        
    case 'maxwell'
        % The input parameters are compliances
        [~,~,~,~,padSize] = makeGeneralizedMaxwellModel(fluidSetting_output,elasticSetting_output,n_terms,timeVec,dt,st);
        [E_storage,E_loss,~,~,lb,ub] = makeHarmonicGeneralizedMaxwellModel(fluidSetting_output,elasticSetting_output,n_terms,timeVec,minTimescale,padSize,forwardFitTimescale);
        
        if strcmp(elasticSetting_output,'y')
            if strcmp(fluidSetting_output,'y')
                tauInds = (3:2:length(ub)-1);
                modulusInds = horzcat(1,(2:2:length(ub)-1));
            else
                tauInds = (3:2:length(ub));
                modulusInds = horzcat(1,(2:2:length(ub)));
            end
        else
            if strcmp(fluidSetting_output,'y')
                tauInds = (2:2:length(ub)-1);
                modulusInds = (1:2:length(ub)-1);
            else
                tauInds = (2:2:length(ub));
                modulusInds = (1:2:length(ub));
            end
        end
        
        % Free up the relaxation times for proper adjustment
        lb(tauInds) = 10^(floor(log10(dt))-1);
        ub(tauInds) = 10^(ceil(log10(timeVec(end)))+1);
        
        % The moduli should still be reasonable, so leave those indices
        % alone
        
        % Make U(s)
        UArms = cell(1,ceil(n_terms + padSize));
        for iii = 1:size(UArms,2)
            if iii == 1 && elasticSetting_input == 'y'
                UArms{iii} = sprintf('c(1)'); % J_g
            elseif iii == size(UArms,2) && fluidSetting_input == 'y'
                UArms{iii} = sprintf('0'); % phi
            else
                if strcmp(elasticSetting_input, 'y')
                    UArms{iii} = sprintf('(c(%d)) ./ (1.0 + ( c(%d).*s ) )',(iii-1)*2,(iii-1)*2+1);
                else
                    UArms{iii} = sprintf('(c(%d)) ./ (1.0 + ( c(%d).*s ) )',(iii)*2-1,(iii)*2);
                end
            end
        end

        % Create the Maxwell storage element series
        tempString1 = '(';
        for iii = 1:size(UArms,2)
            if iii < size(UArms,2)
                tempString1 = horzcat(tempString1, sprintf('(%s)+',UArms{iii}));
            else
                tempString1 = horzcat(tempString1, sprintf('(%s))',UArms{iii}));
            end
        end
        
        UString = sprintf('@(s,c) %s',tempString1);
        Ufunc = str2func(UString);
        
        % Make Q(s)
        QArms = cell(1,ceil(n_terms + padSize));
        for iii = 1:size(QArms,2)
            if iii == 1 && elasticSetting_output == 'y'
                QArms{iii} = sprintf('c(1)'); % E_e
            elseif iii == size(QArms,2) && fluidSetting_output == 'y'
                QArms{iii} = sprintf('0'); % phi
            else
                if strcmp(elasticSetting_output, 'y')
                    QArms{iii} = sprintf('(c(%d).*c(%d).*s) ./ (1.0 + ( c(%d).*s ) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
                else
                    QArms{iii} = sprintf('(c(%d).*c(%d).*s) ./ (1.0 + ( c(%d).*s ) )',(iii)*2-1,(iii)*2,(iii)*2);
                end
            end
        end

        % Create the Maxwell storage element series
        tempString1 = '(';
        for iii = 1:size(QArms,2)
            if iii < size(QArms,2)
                tempString1 = horzcat(tempString1, sprintf('(%s)+',QArms{iii}));
            else
                tempString1 = horzcat(tempString1, sprintf('(%s))',QArms{iii}));
            end
        end
        
        QString = sprintf('@(s,c) %s',tempString1);
        Qfunc = str2func(QString);
        
        % Create Y data and Models
        y_data = Ufunc(slog,inputParams);
        outFunc = @(s,c) 1./(Qfunc(s,c));
        
%         figure
%         plot((slog'),(y_data'),'*')
%         hold on
%         plot((slog'),(objFunc(slog,beta0)),'-')
%         set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         xlabel('s [1/s]')
%         ylabel('U [1/Pa]')
%         hold off
        
        eflag = 0;
        
        beta0 = NaN(size(ub));
%         beta0(tauInds) = 10.^(log10(lb(tauInds))+rand(size(tauInds)).*(log10(ub(tauInds))-log10(lb(tauInds))));
%         beta0(modulusInds) = 10.^(log10(lb(modulusInds))+rand(size(modulusInds)).*(log10(ub(modulusInds))-log10(lb(modulusInds))));
%         if modulusInds(1) == 1
%             beta0(1) = 1./(inputParams(1));
%         end
        beta0(tauInds) = inputParams(tauInds);
        beta0(modulusInds) = 10.^(log10(1./(10.^log10(abs(max(y_data,[],'all')))))+rand(size(modulusInds)).*(log10(1./(10.^log10(abs(min(y_data,[],'all')))))-log10(1./(10.^log10(abs(max(y_data,[],'all')))))));
        if modulusInds(1) == 1
            beta0(1) = 1./(inputParams(1));
        end
        
        runNum = 0;
        bestRes = Inf;
        bestParams = NaN(size(ub));
        
        if strcmp(loopType,'while')
            while any(eflag == [0 2 3])
                runNum = runNum + 1;

                objFunc = @(c,s) vertcat(real(outFunc(s,c)),imag(outFunc(s,c)));
                [beta,resnorm,~,eflag,~,~,~] = ...
                    lsqcurvefit(objFunc,beta0,slog,vertcat(real(y_data),imag(y_data)),lb,[],lsqoptions);

    %             costFunc = @(c) sum(abs((outFunc(slog,c)-y_data)./y_data),'all');
    %             [beta,~,eflag,~] = ...
    %                 fmincon(costFunc,beta0,[],[],[],[],lb,[],[],opts);
    %             resnorm = sum(abs((outFunc(slog,beta)-y_data)./y_data),'all');

                % Show the performance
    %             figure
    %             plot(log10(slog),log10(y_data),'*')
    %             hold on
    %             plot(log10(slog),log10(outFunc(slog,beta0)),'b')
    %             plot(log10(slog),log10(outFunc(slog,beta)),'r')
    %             hold off

                beta0 = NaN(size(ub));
    %             beta0(tauInds) = 10.^(log10(lb(tauInds))+rand(size(tauInds)).*(log10(ub(tauInds))-log10(lb(tauInds))));
    %             beta0(modulusInds) = 10.^(log10(lb(modulusInds))+rand(size(modulusInds)).*(log10(ub(modulusInds))-log10(lb(modulusInds))));
                beta0(tauInds) = inputParams(tauInds);
                beta0(modulusInds) = 10.^(log10(1./(10.^log10(abs(max(y_data,[],'all')))))+rand(size(modulusInds)).*(log10(1./(10.^log10(abs(min(y_data,[],'all')))))-log10(1./(10.^log10(abs(max(y_data,[],'all')))))));
                if modulusInds(1) == 1
                    beta0(1) = 1./(inputParams(1));
                end

                if resnorm < bestRes
                    bestRes = resnorm;
                    bestParams = beta;
                    if showPlots
                        fprintf('New Best Params!\n');
                    end
                end

                if showPlots
                    fprintf('Iteration %d; Exit Flag %d\n',runNum,eflag);
                end

                if runNum >= maxNum
                    if showPlots
                        fprintf('\nMaximum No. of Translation Attempts Reached...The New Predicted Parameters May Be Inaccurate.\n');
                    end
                    break;
                end
            end
        else
            
            if isempty(gcp('nocreate'))
                %parpool(N_workers)
                parpool(8,'IdleTimeout', Inf)
            end
            
            beta_dist = zeros(length(ub),maxNum);
            beta0_dist = zeros(size(beta_dist));
            resnorm_dist = zeros(1,maxNum);
            objFunc = @(c,s) vertcat(real(outFunc(s,c)),imag(outFunc(s,c)));
            
            progressString = sprintf('Converting GKV to GM\nInvestigating S-Plane Functions\nParallel Search Running...');
            hbar = parfor_progressbar(maxNum,progressString);
            warning('off');
            
            parfor i = 1:maxNum
                
                beta0 = NaN(size(ub));
                beta0(tauInds) = 10.^(log10(lb(tauInds))+rand(size(tauInds)).*(log10(ub(tauInds))-log10(lb(tauInds))));
                beta0(modulusInds) = 10.^(log10(lb(modulusInds))+rand(size(modulusInds)).*(log10(ub(modulusInds))-log10(lb(modulusInds))));
%                 beta0(tauInds) = inputParams(tauInds);
%                 beta0(modulusInds) = 10.^(log10(1./(10.^log10(abs(max(y_data,[],'all')))))+rand(size(modulusInds)).*(log10(1./(10.^log10(abs(min(y_data,[],'all')))))-log10(1./(10.^log10(abs(max(y_data,[],'all')))))));
%                e if modulusInds(1) == 1
%                     beta0(1) = 1./(inputParams(1));
%                 end
                
                [beta,resnorm,~,~,~,~,~] = ...
                    lsqcurvefit(objFunc,beta0,slog,vertcat(real(y_data),imag(y_data)),lb,[],lsqoptions);
                
                % Show the performance
%                 figure
%                 plot(log10(slog),log10(y_data),'*')
%                 hold on
% %                 plot(log10(slog),log10(outFunc(slog,beta0)),'b')
%                 plot(log10(slog),log10(outFunc(slog,beta)),'r')
%                 hold off
                
                beta_dist(:,i) = beta;
                beta0_dist(:,i) = beta0;
                resnorm_dist(i) = resnorm;
                hbar.iterate(1) % Increase progressbar by 1 iteration
            end
            
            close(hbar)
            warning('on');
            
            [~,idx] = min(resnorm_dist(resnorm_dist > 0));
            bestBeta0 = beta0_dist(:,idx);
            bestParams = beta_dist(:,idx);
            
        end
        
        beta = bestParams;
        
        if showPlots
            % Show the performance
            fittingResultTest = figure('position',[500 20 1000 1000]);
            plot(log10(slog),log10(y_data),'*')
            hold on
            if strcmp(loopType,'parfor')
%                 plot(log10(slog),log10(outFunc(slog,bestBeta0)),'b')
            end
            plot(log10(slog),log10(outFunc(slog,beta)),'r')
            title('Fitting Result: Conversion to GM', 'FontSize', 16)
            set(findall(gcf,'-property','FontSize'),'FontSize',34)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',5)
            xlabel('$s, log(1/s)$', 'FontSize', 40, 'Interpreter', 'Latex')
            ylabel('$U(s), log(1/Pa)$', 'FontSize', 40, 'Interpreter', 'Latex')
            hold off
            
            % First, generate a frequency array in log scale
            de0 = 2.0.*pi.*(1./timeVec(end));
            maxi = 2.0.*pi.*(1./dt);
            omega = log_tw(de0,maxi);
            
            Jstorage = J_storage_advanced(omega,inputParams,elasticSetting_input,fluidSetting_input);
            Jloss = J_loss_advanced(omega,inputParams,elasticSetting_input,fluidSetting_input);
            Jabs = sqrt(Jstorage.^2+Jloss.^2);
            
            if exist('storagePlotTest','var')
                clearvars storagePlot
            end
            if exist('lossModulusPlotTest','var')
                clearvars lossModulusPlot
            end
            storagePlotTest = figure('position',[500 30 1000 1000]);
            hold on
            grid on
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            title('Storage Modulus: Conversion to GM', 'FontSize', 16)
            set(findall(gcf,'-property','FontSize'),'FontSize',34)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',5)
            xlabel('$\omega, \,rad/s$', 'FontSize', 40, 'Interpreter', 'Latex')
            ylabel('$E''(\omega)$', 'FontSize', 40, 'Interpreter', 'Latex')
            % xlim([1./10^ceil(log10(dataStruct(numFiles+kk).t_r(end))) 1./10^ceil(log10(dataStruct(numFiles+kk).dt))])
            hold off

            lossModulusPlotTest = figure('position',[500 40 1000 1000]);
            hold on
            grid on
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            title('Loss Modulus: Conversion to GM', 'FontSize', 16)
            set(findall(gcf,'-property','FontSize'),'FontSize',34)
            set(gca,'TickLength',[0.02 0.02],'LineWidth',5)
            xlabel('$\omega, \,rad/s$', 'FontSize', 40, 'Interpreter', 'Latex')
            ylabel('$E''''(\omega)$', 'FontSize', 40, 'Interpreter', 'Latex')
            % xlim([1./10^ceil(log10(dataStruct(numFiles+kk).t_r(end))) 1./10^ceil(log10(dataStruct(numFiles+kk).dt))])
            hold off
            
            % Calculate E_storage
            figure(storagePlotTest);
            hold on
            plot(omega, E_storage(omega,beta), '-r', 'LineWidth', 7)
            scatter(omega, Jstorage./(Jabs.^2),markerSize,'o','Linewidth',linewidthSize/2)
            hold off

            % Calculate E_loss
            figure(lossModulusPlotTest);
            hold on
            plot(omega, E_loss(omega,beta), '-r', 'LineWidth', 7)
            scatter(omega, Jloss./(Jabs.^2),markerSize,'o','Linewidth',linewidthSize/2)
            hold off
        end
        
    otherwise
        error('convertParams: There is no feature built to convert parameters to your desired model');

end

end

