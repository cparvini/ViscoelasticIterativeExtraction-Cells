function [U_func,F_conv,F_conv_wrapper,lb,ub,subref,selector] = makeGeneralizedVoigtModel(elasticSetting,fluidSetting,n_terms,timeVec,dt,st,minTimescale,forwardFitTimescale)
%makeGeneralizedVoigtModel Create a Generalized Voigt function of variable
%length
%   This function takes a variety of settings and a number of terms, and
%   correspondingly creates a Generalized Voigt viscoelastic series model
%   using Matlab anonymous functions. U_func is the retardance spectrum of
%   the model, F_conv is the convolution of the force and the retardance.
%   In addition, lb and ub are the lower and upper bounds respectively that
%   are used during nonlinear least-squares fitting to appropriately limit
%   the parameters while they are varied by the solver.
    
    % Performing a convolution results in data that is N+M+1 in length. In
    % order to continue, we are selecting the first portion of this result,
    % which is the relevant portion for our analysis.
    subref = @(A) A(1:size(timeVec,1),1);
    
    % In the special case of an advanced log approach, we have to specify
    % only one of the outputs is used and deal it to the resulting variable
    % instead.
    selector = @(x) deal(x(1,:));
    
    % Limits for the parameters
    lbLimit = 1e-12;
    ubLimit = 10;
    convType = 'full';

    % Function Shape for Multiple Voigt with No Load Assumption
    % Create list of compliances
    if fluidSetting == 'y' && elasticSetting == 'y'
        compArms = cell(1,(n_terms+2));
        lb = zeros(1,n_terms*2+2);
        ub = ones(1,n_terms*2+2);
        initParams = zeros(1,n_terms*2+2);
    elseif fluidSetting == 'n' && elasticSetting == 'n'
        compArms = cell(1,(n_terms));
        lb = zeros(1,n_terms*2);
        ub = ones(1,n_terms*2);
        initParams = zeros(1,n_terms*2);
    else
        compArms = cell(1,(n_terms+1));
        lb = zeros(1,n_terms*2+1);
        ub = ones(1,n_terms*2+1);
        initParams = zeros(1,n_terms*2+1);
    end
    
    lb = lb + lbLimit;
    ub = ub * ubLimit;

    for iii = 1:size(compArms,2)
        if iii == 1 && elasticSetting == 'y'
            compArms{iii} = sprintf('0');
        elseif iii == size(compArms,2) && fluidSetting == 'y'
            compArms{iii} = sprintf('c(%d)',(iii-1)*2); % Changed to match Lopez et al. code (Github, Lib_rheology.py)
        else
            if strcmp(elasticSetting, 'y')
                compArms{iii} = sprintf('(c(%d)./c(%d)).*exp( (-inputs(:,1))./c(%d) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
            else
                compArms{iii} = sprintf('(c(%d)./c(%d)).*exp( (-inputs(:,1))./c(%d) )',(iii)*2-1,(iii)*2,(iii)*2);
            end
        end
    end

    % Create the voigt element series
    tempString1 = '(';
    for iii = 1:size(compArms,2)
        if iii < size(compArms,2)
            tempString1 = horzcat(tempString1, sprintf('(%s)+',compArms{iii}));
        else
            tempString1 = horzcat(tempString1, sprintf('(%s))',compArms{iii}));
        end
    end

    if n_terms >= 1
        if strcmp(convType,'full')
            compString = '@(c,inputs) (convnfft(( ';
            compString = horzcat(compString, sprintf('%s',tempString1) );
            compString = horzcat(compString, sprintf('),inputs(:,2),''full''))') );
            retardanceString = sprintf('@(c,inputs) %s',tempString1);
            U_func_temp = str2func(retardanceString);
            F_conv_temp = str2func(compString);
            if strcmp(elasticSetting, 'y')
                U_func = @(c,inputs) c(1)+U_func_temp(c,inputs);
                F_conv = @(c,inputs) c(1).*inputs(:,2)+subref(F_conv_temp(c,inputs))*dt;
            else
                U_func = @(c,inputs) U_func_temp(c,inputs);
                F_conv = @(c,inputs) subref(F_conv_temp(c,inputs))*dt;
            end
        else
            compString = '@(c,inputs) convnfft(( ';
            compString = horzcat(compString, sprintf('%s',tempString1) );
            compString = horzcat(compString, sprintf('),inputs(:,2),''same'')') );
            retardanceString = sprintf('@(c,inputs) %s',tempString1);
            U_func_temp = str2func(retardanceString);
            F_conv_temp = str2func(compString);
            if strcmp(elasticSetting, 'y')
                U_func = @(c,inputs) c(1)+U_func_temp(c,inputs);
                F_conv = @(c,inputs) c(1).*inputs(:,2)+(F_conv_temp(c,inputs))*dt;
            else
                U_func = @(c,inputs) U_func_temp(c,inputs);
                F_conv = @(c,inputs) (F_conv_temp(c,inputs))*dt;
            end
        end
    else
        compString = '@(c,inputs) zeros(size(inputs(:,2)))';
        F_conv = str2func(compString);
        if strcmp(elasticSetting, 'y')
            retardanceString = sprintf('@(c,inputs) c(1).*ones(size(inputs(:,2)))');
            U_func = str2func(retardanceString);
        else
            retardanceString = '@(c,inputs) zeros(size(inputs(:,2)))';
            U_func = str2func(retardanceString);
        end
    end

    F_conv_wrapper = @(c,inputs) log_scale(...
        (F_conv(c,inputs)),...
        inputs(:,1),dt,st);
        
    timescaleArray = [];
    for iii = 1:size(compArms,2)
    
        if iii == 1 && elasticSetting == 'y'
            % Initial guess for glassy compliance
            initParams(1) = 1e-6;
            continue;
        end

        if iii == size(compArms,2) && fluidSetting == 'y'
            % Initial guess, and bounds, for fluidity
            initParams(end) = 10^0;
            lb(end) = 0;
            ub(end) = ubLimit;
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

            initParams(iii1) = 1e-6;
            initParams(iii2) = (10^(iii-1))*minTimescale;
        else
            topTimescale = ceil(log10(max(timeVec)))+1;
            
            initParams(iii1) = 1e-6;
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
            scaledLB = logspace(topTimescale-1,floor(log10(minTimescale)),n_terms)...
                * ( minTimescale / 10^floor(log10(minTimescale)) );
            lb(timescaleArray) = scaledLB;
            initParams(timescaleArray) = lb(timescaleArray).*10;
            ub(timescaleArray) = lb(timescaleArray).*100;
        end
    end

end

