function [J_biprime] = J_loss_advanced(omega, params, elasticSetting, fluidSetting)
%J_loss  This function gives an array of loss compliance on omega.
%   This function calculates the loss compliance given an input
%   frequency vector, omega, and the model parameters contained within
%   params.

    J_biprime = zeros(size(omega));
    
    if ischar(elasticSetting)
       if strcmp(elasticSetting,'y')
            elasticSetting = 1;
       else
            elasticSetting = 0; 
       end
    end
    
    if ischar(fluidSetting)
       if strcmp(fluidSetting,'y')
            fluidSetting = 1;
       else
            fluidSetting = 0; 
       end
    end
    
    if elasticSetting && ~fluidSetting
        J = params(2:2:end);
        tau = params(3:2:end);
        phi = 0;

        for jj = 1:length(omega)
            if length(J) > 1
                J_biprime(jj) = sum( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            else
                J_biprime(jj) = ( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            end
        end

    elseif ~elasticSetting && fluidSetting
        J = params(1:2:(end-1));
        tau = params(2:2:(end-1));
        phi = params(end);

        for jj = 1:length(omega)
            if length(J) > 1
                J_biprime(jj) = sum( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            else
                J_biprime(jj) = ( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            end
        end

    elseif elasticSetting && fluidSetting
        J = params(2:2:(end-1));
        tau = params(3:2:(end-1));
        phi = params(end);

        for jj = 1:length(omega)
            if length(J) > 1
                J_biprime(jj) = sum( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            else
                J_biprime(jj) = ( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            end
        end

    else
        J = params(1:2:end);
        tau = params(2:2:end);
        phi = 0;

        for jj = 1:length(omega)
            if length(J) > 1
                J_biprime(jj) = sum( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            else
                J_biprime(jj) = ( (J .* omega(jj) .* tau) ./(1.0 + ((omega(jj).^2).*(tau.^2)) ) ) + phi / omega(jj);
            end
        end  
            
    end

end

