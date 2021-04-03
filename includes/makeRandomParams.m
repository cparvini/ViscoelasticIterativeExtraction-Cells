function [beta0,tauInds,modulusInds] = makeRandomParams(oldParams,ub_randLimit,lb_randLimit,elasticSetting,fluidSetting,varargin)
%makeRandomParams Create Random Parameters Within Bounds
%   This function takes in the relevant elastic and fluid settings for the
%   viscoelastic model, including the limits within which random guesses
%   should be generated, and provides both a random starting point the
%   indices that correspond to characteristic times (tauInds) or
%   stiffnesses (modulusInds).

    % Initialize
    newInds = (false(size(ub_randLimit)));
    enforceGridGuess = 0;
    elasticInd = 1;
    
    % Checking varargin structure
    if ~isempty(varargin)
        if length(varargin) > 1
            for i = 1:length(varargin)
               switch i
                   case 1
                       newInds = (cell2mat(varargin{i}));
                   case 2
                       enforceGridGuess = logical(cell2mat(varargin{i}));
                   case 3
                       ub_loop = (cell2mat(varargin{i}));
                   case 4
                       lb_loop = (cell2mat(varargin{i}));
                   case 5
                       loopind = (cell2mat(varargin{i}));
                   otherwise
                       error('Too many arguments passed to makeRandomParams via varargin!')
               end
            end
        else
            newInds = (cell2mat(varargin));
        end
    end
    
    beta0 = oldParams;

    beta0_temp = zeros(size(ub_randLimit));
    if any(newInds) && ~enforceGridGuess
        
        if strcmp(elasticSetting,'y')
            if strcmp(fluidSetting,'y')
                beta0_temp(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                beta0_temp(2:2:end-1) = 10.^(((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1)))+lb_randLimit(2:2:end-1)));
                beta0_temp(3:2:end-1) = ((ub_randLimit(3:2:end-1)-lb_randLimit(3:2:end-1)).*rand(size(ub_randLimit(3:2:end-1))) + lb_randLimit(3:2:end-1));
                beta0_temp(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                tauInds = (3:2:length(beta0_temp)-1);
                modulusInds = horzcat(elasticInd,(2:2:length(beta0_temp)-1));
            else
                beta0_temp(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                beta0_temp(2:2:end) = 10.^(((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end)))+lb_randLimit(2:2:end)));
                beta0_temp(3:2:end) = ((ub_randLimit(3:2:end)-lb_randLimit(3:2:end)).*rand(size(ub_randLimit(3:2:end))) + lb_randLimit(3:2:end));
                tauInds = (3:2:length(beta0_temp));
                modulusInds = horzcat(elasticInd,(2:2:length(beta0_temp)));
            end
        else
            if strcmp(fluidSetting,'y')
                beta0_temp(1:2:end-1) = 10.^(((ub_randLimit(1:2:end-1)-lb_randLimit(1:2:end-1)).*rand(size(ub_randLimit(1:2:end-1)))+lb_randLimit(1:2:end-1)));
                beta0_temp(2:2:end-1) = ((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1))) + lb_randLimit(2:2:end-1));
                beta0_temp(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                tauInds = (2:2:length(beta0_temp)-1);
                modulusInds = (1:2:length(beta0_temp)-1);
            else
                beta0_temp(1:2:end) = 10.^(((ub_randLimit(1:2:end)-lb_randLimit(1:2:end)).*rand(size(ub_randLimit(1:2:end)))+lb_randLimit(1:2:end)));
                beta0_temp(2:2:end) = ((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end))) + lb_randLimit(2:2:end));
                tauInds = (2:2:length(beta0_temp));
                modulusInds = (1:2:length(beta0_temp));
            end
        end

    elseif any(newInds) && enforceGridGuess

        if strcmp(elasticSetting,'y')
            if strcmp(fluidSetting,'y')
                tauInds = (3:2:length(ub_loop)-1);
                modulusInds = horzcat(elasticInd,(2:2:length(ub_loop)-1));
            else
                tauInds = (3:2:length(ub_loop));
                modulusInds = horzcat(elasticInd,(2:2:length(ub_loop)));
            end
        else
            if strcmp(fluidSetting,'y')
                tauInds = (2:2:length(ub_loop)-1);
                modulusInds = (1:2:length(ub_loop)-1);
            else
                tauInds = (2:2:length(ub_loop));
                modulusInds = (1:2:length(ub_loop));
            end
        end

        newTau = tauInds(ismember(tauInds,find(newInds)));
        newMod = modulusInds(ismember(modulusInds,find(newInds)));

        tauStep = floor(sqrt(loopLim));
        modStep = floor((loopLim)^(1./(length(ub_loop(newMod)))))/(length(newTau)*tauStep);
        [betaGridJ,betaGridTau] = meshgrid(logspace(log10(lb_loop(newMod)),log10(ub_loop(newMod)),modStep),...
            logspace(log10(lb_loop(newTau)),log10(ub_loop(newTau)),tauStep));

        if loopind <= tauStep*modStep
            beta0_temp(newInds) = [betaGridJ(loopind);betaGridTau(loopind)];
        else
            beta0_temp(newInds) = [betaGridJ(end);betaGridTau(end)];
        end

    else

        if strcmp(elasticSetting,'y')
            if strcmp(fluidSetting,'y')
                beta0_temp(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                beta0_temp(2:2:end-1) = 10.^(((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1)))+lb_randLimit(2:2:end-1)));
                beta0_temp(3:2:end-1) = ((ub_randLimit(3:2:end-1)-lb_randLimit(3:2:end-1)).*rand(size(ub_randLimit(3:2:end-1))) + lb_randLimit(3:2:end-1));
                beta0_temp(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                tauInds = (3:2:length(beta0_temp)-1);
                modulusInds = horzcat(elasticInd,(2:2:length(beta0_temp)-1));
            else
                beta0_temp(1) = 10^(((ub_randLimit(1)-lb_randLimit(1))*rand()+lb_randLimit(1)));
                beta0_temp(2:2:end) = 10.^(((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end)))+lb_randLimit(2:2:end)));
                beta0_temp(3:2:end) = ((ub_randLimit(3:2:end)-lb_randLimit(3:2:end)).*rand(size(ub_randLimit(3:2:end))) + lb_randLimit(3:2:end));
                tauInds = (3:2:length(beta0_temp));
                modulusInds = horzcat(elasticInd,(2:2:length(beta0_temp)));
            end
        else
            if strcmp(fluidSetting,'y')
                beta0_temp(1:2:end-1) = 10.^(((ub_randLimit(1:2:end-1)-lb_randLimit(1:2:end-1)).*rand(size(ub_randLimit(1:2:end-1)))+lb_randLimit(1:2:end-1)));
                beta0_temp(2:2:end-1) = ((ub_randLimit(2:2:end-1)-lb_randLimit(2:2:end-1)).*rand(size(ub_randLimit(2:2:end-1))) + lb_randLimit(2:2:end-1));
                beta0_temp(end) = 10^(floor(log10(ub_fluidityRandLimit_init))*rand());
                tauInds = (2:2:length(beta0_temp)-1);
                modulusInds = (1:2:length(beta0_temp)-1);
            else
                beta0_temp(1:2:end) = 10.^(((ub_randLimit(1:2:end)-lb_randLimit(1:2:end)).*rand(size(ub_randLimit(1:2:end)))+lb_randLimit(1:2:end)));
                beta0_temp(2:2:end) = ((ub_randLimit(2:2:end)-lb_randLimit(2:2:end)).*rand(size(ub_randLimit(2:2:end))) + lb_randLimit(2:2:end));
                tauInds = (2:2:length(beta0_temp));
                modulusInds = (1:2:length(beta0_temp));
            end
        end

    end
    
    if any(newInds)
        beta0(newInds) = beta0_temp(newInds);
    else
        beta0 = beta0_temp;
    end
    
end

