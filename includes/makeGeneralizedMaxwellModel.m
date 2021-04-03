function [G_func,Q_func,h_conv_maxwell,h_conv_maxwell_wrapper,padSize] = makeGeneralizedMaxwellModel(fluidSetting,elasticSetting,n_terms,timeVec,dt,st)
%makeGeneralizedMaxwellModel Create the Generalized Maxwell Model
%   This function creates the force convolution form of the Generalized
%   Maxwell Model in terms of modulus (the inverse of compliance). This
%   function will be used to compare the performance of the parameter sets
%   that are fit to the harmonic quantities from the Generalized Voigt
%   to the actual data. F_conv_maxwell gives the time domain form of the
%   predicted force convolution, and F_conv_maxwell_wrapper gives the log
%   scaled form of the predicted force convolution.

    % Performing a convolution results in data that is N+M+1 in length. In
    % order to continue, we are selecting the first portion of this result,
    % which is the relevant portion for our analysis.
    subref = @(A) A(1:size(timeVec,1),1);
    convType = 'full';
    if (~strcmp(convType,'full') && ~strcmp(convType,'same')) error('The convType inside makeGeneralizedMaxwellModel() is not set properly. Please choose either "full" or "same".'); end
    modType = 'elastic'; % Choose 'elastic' or 'glassy' for model construction
    if (~strcmp(modType,'elastic') && ~strcmp(modType,'glassy')) error('The modType inside makeGeneralizedMaxwellModel() is not set properly. Please choose either "glassy" or "elastic".'); end
    relaxForm = 'analytical'; % Choose 'analytical' or 'discrete' for generation of Q(t)
    if (~strcmp(relaxForm,'analytical') && ~strcmp(relaxForm,'discrete')) error('The relaxForm inside makeGeneralizedMaxwellModel() is not set properly. Please choose either "analytical" or "discrete".'); end
    
    % Create list of Moduli
    if fluidSetting == 'n' && elasticSetting == 'n'
        stiffArmsMaxwell = cell(1,(n_terms));
        padSize = 0;
    else
        stiffArmsMaxwell = cell(1,(n_terms+1));
        padSize = 1;
    end
    
    relaxArmsMaxwell = cell(size(stiffArmsMaxwell));

    % Create the Maxwell relaxance model
    % Create list of stiffnesses 
    for iii = 1:size(stiffArmsMaxwell,2)
        if iii == 1 && elasticSetting == 'y' && fluidSetting == 'n'
            stiffArmsMaxwell{iii} = sprintf('0'); % Included down below
            relaxArmsMaxwell{iii} = sprintf('0'); % Included down below
        elseif iii == 1 && elasticSetting == 'y' && fluidSetting == 'y'
            stiffArmsMaxwell{iii} = sprintf('c(%d).*exp( (-inputs(:,1))./c(end) )',1);
            relaxArmsMaxwell{iii} = sprintf('(-c(%d)./c(end)).*exp( (-inputs(:,1))./c(end) )',1);
        else
            if strcmp(elasticSetting, 'y')
                if strcmp(modType,'elastic')
                    stiffArmsMaxwell{iii} = sprintf('c(%d).*exp( (-inputs(:,1))./c(%d) )',(iii-1)*2,(iii-1)*2+1);
                    relaxArmsMaxwell{iii} = sprintf('-(c(%d)./c(%d)).*exp( (-inputs(:,1))./c(%d) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
                else            
                    stiffArmsMaxwell{iii} = sprintf('(-c(%d)).*(1 - exp((-inputs(:,1))./c(%d)) )',(iii-1)*2,(iii-1)*2+1);
                    relaxArmsMaxwell{iii} = sprintf('-(c(%d)./c(%d)).*(exp((-inputs(:,1))./c(%d)) )',(iii-1)*2,(iii-1)*2+1,(iii-1)*2+1);
                end
            else
                if strcmp(modType,'elastic')
                    stiffArmsMaxwell{iii} = sprintf('c(%d).*exp( (-inputs(:,1))./c(%d) )',(iii)*2-1,(iii)*2);
                    relaxArmsMaxwell{iii} = sprintf('-(c(%d)./c(%d)).*exp( (-inputs(:,1))./c(%d) )',(iii)*2-1,(iii)*2,(iii)*2);
                else
                    stiffArmsMaxwell{iii} = sprintf('(-c(%d)).*(1 - exp((-inputs(:,1))./c(%d)) )',(iii)*2-1,(iii)*2);
                    relaxArmsMaxwell{iii} = sprintf('-(c(%d)./c(%d)).*(exp((-inputs(:,1))./c(%d)) )',(iii)*2-1,(iii)*2,(iii)*2);
                end
            end
        end
    end
    
    [~,modulusInds] = getParamIndices(ones(2*n_terms+strcmp(elasticSetting, 'y')+strcmp(fluidSetting, 'y'),1),elasticSetting,fluidSetting);

    % Create the relaxance series
    tempString1 = '(';
    tempString2 = '(';
    for iii = 1:(n_terms+padSize)
%         if iii == 1 && elasticSetting == 'y'
%             continue;
%         end
        if iii < n_terms+padSize
            tempString1 = horzcat(tempString1, sprintf('(%s)+',stiffArmsMaxwell{iii}));
            tempString2 = horzcat(tempString2, sprintf('(%s)+',relaxArmsMaxwell{iii}));
        else
            tempString1 = horzcat(tempString1, sprintf('(%s))',stiffArmsMaxwell{iii}));
            tempString2 = horzcat(tempString2, sprintf('(%s))',relaxArmsMaxwell{iii}));
        end
    end

    if n_terms >= 1
        if strcmp(convType,'full')
            stiffString = sprintf('@(c,inputs) (%s)',tempString1);
            relaxString = sprintf('@(c,inputs) (%s)',tempString2);
            G_func_temp = str2func(stiffString);
            Q_func_temp = str2func(relaxString);
            if strcmp(elasticSetting, 'y') && fluidSetting == 'n'
                if strcmp(modType,'glassy')
                    if strcmp(relaxForm,'discrete')
                        h_conv_maxwell = @(c,inputs) c(1).*(inputs(:,2))-subref(convnfft(gradient(G_func_temp(c,inputs),dt), inputs(:,2),'full'))*dt;
                    else
                        h_conv_maxwell = @(c,inputs) c(1).*(inputs(:,2))+subref(convnfft(Q_func_temp(c,inputs), inputs(:,2),'full'))*dt;
                    end
                else
                    if strcmp(relaxForm,'discrete')
                        h_conv_maxwell = @(c,inputs) sum(c(modulusInds)).*(inputs(:,2))-subref(convnfft(gradient(G_func_temp(c,inputs),dt), inputs(:,2),'full'))*dt;
                    else
                        h_conv_maxwell = @(c,inputs) sum(c(modulusInds)).*(inputs(:,2))+subref(convnfft(Q_func_temp(c,inputs), inputs(:,2),'full'))*dt;
                    end
                end
            else
                if strcmp(relaxForm,'discrete')
                    h_conv_maxwell = @(c,inputs) -subref(convnfft(gradient(G_func_temp(c,inputs),dt), inputs(:,2),'full'))*dt;
                else
                    h_conv_maxwell = @(c,inputs) subref(convnfft(Q_func_temp(c,inputs), inputs(:,2),'full'))*dt;
                end
            end
        else
            stiffString = sprintf('@(c,inputs) (%s)',tempString1);
            relaxString = sprintf('@(c,inputs) (%s)',tempString2);
            G_func_temp = str2func(stiffString);
            Q_func_temp = str2func(relaxString);
            if strcmp(elasticSetting, 'y') && fluidSetting == 'n'
                if strcmp(modType,'glassy')
                    if strcmp(relaxForm,'discrete')
                        h_conv_maxwell = @(c,inputs) c(1).*(inputs(:,2))-(convnfft(gradient(G_func_temp(c,inputs),dt), inputs(:,2),'same'))*dt;
                    else
                        h_conv_maxwell = @(c,inputs) c(1).*(inputs(:,2))+(convnfft(Q_func_temp(c,inputs), inputs(:,2),'same'))*dt;
                    end
                else
                    if strcmp(relaxForm,'discrete')
                        h_conv_maxwell = @(c,inputs) sum(c(modulusInds)).*(inputs(:,2))-(convnfft(gradient(G_func_temp(c,inputs),dt), inputs(:,2),'same'))*dt;
                    else
                        h_conv_maxwell = @(c,inputs) sum(c(modulusInds)).*(inputs(:,2))+(convnfft(Q_func_temp(c,inputs), inputs(:,2),'same'))*dt;
                    end
                end
            else
                if strcmp(relaxForm,'discrete')
                    h_conv_maxwell = @(c,inputs) -(convnfft(gradient(G_func_temp(c,inputs),dt), inputs(:,2),'same'))*dt;
                else
                    h_conv_maxwell = @(c,inputs) (convnfft(Q_func_temp(c,inputs), inputs(:,2),'same'))*dt;
                end
            end
        end
    else
        relaxString = '@(c,inputs) zeros(size(inputs(:,2)))';
        Q_func_temp = str2func(relaxString);
        if strcmp(elasticSetting, 'y') && strcmp(fluidSetting, 'n')
            stiffString = '@(c,inputs) c(1).*ones(size(inputs(:,2)))';
            h_conv_maxwell = str2func(stiffString);
        else
            stiffString = '@(c,inputs) zeros(size(inputs(:,2)))';
            h_conv_maxwell = str2func(stiffString);
        end
    end
    
    G_func = str2func(stiffString);
    if strcmp(relaxForm,'discrete')
        if strcmp(elasticSetting, 'y') && strcmp(fluidSetting, 'n')
            Q_func = @(c,inputs) c(1)-gradient(G_func(c,inputs),dt);
        else
            Q_func = @(c,inputs) -gradient(G_func(c,inputs),dt);
        end
    else
        if strcmp(elasticSetting, 'y') && strcmp(fluidSetting, 'n')
            Q_func = @(c,inputs) c(1)+Q_func_temp(c,inputs);
        else
            Q_func = @(c,inputs) Q_func_temp(c,inputs);
        end
    end
    
    % Write final log function
    h_conv_maxwell_wrapper = @(c,inputs) log_scale(...
        (h_conv_maxwell(c,inputs)),...
        inputs(:,1),dt,st);
    
end

