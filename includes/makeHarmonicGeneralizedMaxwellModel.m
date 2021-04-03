function [E_storage,E_loss,lossAngle,harmonic_wrapper,lb,ub,varargout] = makeHarmonicGeneralizedMaxwellModel(fluidSetting,elasticSetting,n_loss_fit,timeVec,minTimescale,padSize,forwardFitTimescale,varargin)
%makeHarmonicGeneralizedMaxwellModel Create the Generalized Maxwell Model
%in terms of Harmonic Quantities
%   This function creates the harmonic distributions that are characterized
%   by a series of moduli and characteristic times. The anonymous functions
%   created here, E_storage, E_loss, and LossAngle, are all functions of
%   frequency and a parameter set, which is found using the
%   harmonic_wrapper function and a nonlinear least squares fitting process
%   that takes place after this function has been executed. The
%   harmonic_wrapper is written in such a manner that it can be fed
%   directly to an lsqcurvefit call with the lower and upper bounds (lb,
%   ub). It is compared with the inverse of the Loss Compliance and Storage
%   Compliance, which are found after the Generalized Voigt fitting
%   process. The inverse of both compliances are stacked, and compared to
%   the Storage and Loss Modului calculated with E_storage and E_loss such
%   that they can both be simultaneously fit using lsqcurvefit. The result
%   should be a set of parameters (moduli and characteristic times) that
%   match the harmonic quantities calculated from the Generalized Voigt
%   model, which should have in turn given a very good estimate of the
%   force convolution. The function "makeGeneralizedMaxwellModel" is used
%   to calculate the force convolution in time, in conjunction with the 
%   parameter set found by fitting harmonic_wrapper to the Storage and Loss
%   Compliance.

    % Checking varargin structure
    if ~isempty(varargin)
        vararg = cell2mat(varargin{1});
        switch(vararg)
            case 0
                normParams = 'none';
            case 1 
                normParams = 'mod';
            case 2
                normParams = 'tau';
            otherwise
                error('Unknown Command Passed to Maxwell Model Generation (varargin term)');
        end
    end
    
    modType = 'elastic'; % Choose 'elastic' or 'glassy' for model construction
    if (~strcmp(modType,'elastic') && ~strcmp(modType,'glassy')) error('The modType inside makeHarmonicGeneralizedMaxwellModel() is not set properly. Please choose either "glassy" or "elastic".'); end

    storageArms = cell(1,ceil(n_loss_fit + padSize));
    lossArms = cell(1,ceil(n_loss_fit + padSize));
    lbLimit = 0.1;
    ubLimit = 12;
    
    lb = zeros(1,n_loss_fit*2+strcmp(elasticSetting,'y')+strcmp(fluidSetting,'y'))+lbLimit;
    ub = ones(1,n_loss_fit*2+strcmp(elasticSetting,'y')+strcmp(fluidSetting,'y'))*(10^(ubLimit));
    initParams = zeros(1,n_loss_fit*2+padSize);

    timescaleArray = [];
    for iii = 1:size(lossArms,2)

        if iii == 1 && elasticSetting == 'y'
            % Initial guess for glassy compliance
            initParams(1) = 1e6;
            continue;
        end

        if iii == size(lossArms,2) && fluidSetting == 'y'
            % Initial guess, and bounds, for fluidity
            initParams(end) = 10^0;
            lb(end) = 0;
            ub(end) = 10^(ubLimit);
            continue;
        end
        
        if strcmp(elasticSetting, 'y')
            iii1 = ((iii-1)*2);
            iii2 = ((iii-1)*2+1);
        else
            iii1 = (iii*2-1);
            iii2 = (iii*2);
        end

        if forwardFitTimescale
            % Retardation Time
            lb(iii2) = (10^(iii-1))*minTimescale/10;
            ub(iii2) = (10^(iii-1))*minTimescale*10;

            initParams(iii1) = 1e6;
            initParams(iii2) = (10^(iii-1))*minTimescale;
        else
            topTimescale = ceil(log10(max(timeVec)))+1;
            
            initParams(iii1) = 1e6;
            initParams(iii2) = (10^(topTimescale)) / (10^(iii-1)) * ( minTimescale / 10^floor(log10(minTimescale)) );
            
            timescaleArray = horzcat(timescaleArray,iii2);
            
            % Retardation Time
            lb(iii2) = (10^(topTimescale)) / (10^(iii)) * ( minTimescale / 10^floor(log10(minTimescale)) );
            ub(iii2) = lb(iii2)*100;
        end
        
    end
    
    % Check if we need to rescale the timescales to fit within the user
    % defined bounds.
    if ~forwardFitTimescale
        if any(lb(timescaleArray) < minTimescale)
            scaledLB = logspace(topTimescale-1,floor(log10(minTimescale)),n_loss_fit)...
                * ( minTimescale / 10^floor(log10(minTimescale)) );
            lb(timescaleArray) = scaledLB;
            initParams(timescaleArray) = lb(timescaleArray).*10;
            ub(timescaleArray) = lb(timescaleArray).*100;
        end
    end
    
    if strcmp(elasticSetting,'y')
        if strcmp(fluidSetting,'y')
            tauInds = (3:2:length(ub)-1);
            modulusInds = horzcat(1,(2:2:length(ub)-1));
        else
            tauInds = (3:2:length(ub));
            modulusInds = horzcat(1,(2:2:length(ub)));
        end
    else
        if strcmp(fluidSetting,'y')
            tauInds = (2:2:length(ub)-1);
            modulusInds = (1:2:length(ub)-1);
        else
            tauInds = (2:2:length(ub));
            modulusInds = (1:2:length(ub));
        end
    end
    
    if exist('normParams','var')
        switch(normParams)
            case 'none'
                varargout = {1};
%                 varargout = {1,1};
            case 'mod'
                varargout = {abs( 10^(((max(log10(ub(modulusInds)))-min(log10(lb(modulusInds))))./2) - ((max(log10(ub(tauInds)))-min(log10(lb(tauInds))))./2)) )};
%                 varargout = {1,1};
            case 'tau'
                varargout = {1};
%                 varargout = {1,1};
        end
    end
    
    for iii = 1:size(storageArms,2)
        if iii == 1 && elasticSetting == 'y'
            storageArms{iii} = sprintf('c(1)'); % Elastic term, Ee or Eg
        elseif iii == size(storageArms,2) && fluidSetting == 'y'
            storageArms{iii} = sprintf('0'); % SS Fluidity doesn't contribute to storage modulus
        else
            if strcmp(elasticSetting, 'y')
                if strcmp(modType,'elastic')
                    storageArms{iii} = sprintf('(c(%d).*(omega.^2).*(c(%d).^2)) ./ (1.0 + ( (omega.^2).*(c(%d).^2) ) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
                else
                    storageArms{iii} = sprintf('(-c(%d)) ./ (1.0 + ( (omega.^2).*(c(%d).^2) ) )',(iii-1)*2,(iii-1)*2+1);
                end
            else
                if strcmp(modType,'elastic')
                    storageArms{iii} = sprintf('(c(%d).*(omega.^2).*(c(%d).^2)) ./ (1.0 + ( (omega.^2).*(c(%d).^2) ) )',(iii)*2-1,(iii)*2,(iii)*2);
                else
                    storageArms{iii} = sprintf('(-c(%d)) ./ (1.0 + ( (omega.^2).*(c(%d).^2) ) )',(iii-1)*2,(iii-1)*2+1);
                end
            end
        end
    end

    % Create the Maxwell storage element series
    tempString1 = '(';
    for iii = 1:size(storageArms,2)
        if iii < size(storageArms,2)
            tempString1 = horzcat(tempString1, sprintf('(%s)+',storageArms{iii}));
        else
            tempString1 = horzcat(tempString1, sprintf('(%s))',storageArms{iii}));
        end
    end

    
    for iii = 1:size(lossArms,2)
        if iii == 1 && elasticSetting == 'y'
            lossArms{iii} = sprintf('0'); % Elastic Term. If SS Fluidity is included, the elastic term will be included in the last cell entry.
        elseif iii == size(lossArms,2) && elasticSetting == 'y' && fluidSetting == 'y'
%             lossArms{iii} = sprintf('c(%d) ./ omega', (iii-1)*2); % phi
            lossArms{iii} = sprintf('(c(%d).*omega.*c(end)) ./ (1.0 + ( (omega.^2).*(c(end).^2) ) )',1); % tau_f
        elseif iii == size(lossArms,2) && elasticSetting == 'n' && fluidSetting == 'y'
%             lossArms{iii} = sprintf('c(%d) ./ omega', (iii-1)*2); % phi
            lossArms{iii} = sprintf('(omega.*c(end)) ./ (1.0 + ( (omega.^2).*(c(end).^2) ) )'); % tau_f
        else
            if strcmp(elasticSetting, 'y')
                lossArms{iii} = sprintf('(c(%d).*omega.*c(%d)) ./ (1.0 + ( (omega.^2).*(c(%d).^2) ) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
            else
                lossArms{iii} = sprintf('(c(%d).*omega.*c(%d)) ./ (1.0 + ( (omega.^2).*(c(%d).^2) ) )',(iii)*2-1,(iii)*2,(iii)*2);
            end
        end
    end

    % Create the Maxwell loss element series
    tempString2 = '(';
    for iii = 1:size(lossArms,2)
        if iii < size(lossArms,2)
            tempString2 = horzcat(tempString2, sprintf('(%s)+',lossArms{iii}));
        else
            tempString2 = horzcat(tempString2, sprintf('(%s))',lossArms{iii}));
        end
    end

    storageString = sprintf('@(omega,c) %s',tempString1);
    lossString = sprintf('@(omega,c) %s',tempString2);
    lossAngleString = sprintf('@(omega,c) atan(%s./%s).*(180./pi)',tempString2,tempString1);

    E_storage = str2func(storageString);
    E_loss = str2func(lossString);
    lossAngle = str2func(lossAngleString);
    harmonic_wrapper = @(c,input) horzcat(E_loss(input,c),E_storage(input,c));

end

