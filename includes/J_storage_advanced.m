function [J_prime] = J_storage_advanced(omega, params, elasticSetting, fluidSetting)
%J_storage  This function gives an array of storage compliance on omega.
%   This function calculates the storage compliance given an input
%   frequency vector, omega, and the model parameters contained within
%   params.
    
    J_prime = zeros(size(omega));
    
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
        Jg = params(1);
        J = params(2:2:end);
        tau = params(3:2:end);

        for jj = 1:length(omega)
            if length(J) > 1
%                 J_prime(jj) = Jg + sum( J ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (Jg+sum(J)) - sum( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
            else
%                 J_prime(jj) = Jg + ( J ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (Jg+J) - ( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
            end
        end

    elseif ~elasticSetting && fluidSetting
        J = params(1:2:(end-1));
        tau = params(2:2:(end-1));

        for jj = 1:length(omega)
            if length(J) > 1
%                 J_prime(jj) = sum( J ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (sum(J)) - sum( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
            else
%                 J_prime(jj) = ( J ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (J) - ( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
            end
        end

    elseif elasticSetting && fluidSetting
        Jg = params(1);
        J = params(2:2:(end-1));
        tau = params(3:2:(end-1));

        for jj = 1:length(omega)
            if length(J) > 1
%                 J_prime(jj) = Jg + sum( J ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (Jg+sum(J)) - sum( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
            else
%                 J_prime(jj) = Jg + ( J ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (Jg+J) - ( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
            end
        end

    else
        J = params(1:2:end);
        tau = params(2:2:end);

        for jj = 1:length(omega)
            if length(J) > 1
%                 J_prime(jj) = sum( J ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (sum(J)) - sum( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + ( (omega(jj).^2).*(tau.^2) ) ) );
            else
%                 J_prime(jj) = ( J ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
                J_prime(jj) = (J) - ( ( J.*(omega(jj).^2).*(tau.^2) ) ./ (1.0 + (  (omega(jj).^2).*(tau.^2) ) ) );
            end
        end   
            
    end


end

