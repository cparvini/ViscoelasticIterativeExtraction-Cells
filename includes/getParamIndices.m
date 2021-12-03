function [tauInds,modulusInds] = getParamIndices(beta,elasticSetting,fluidSetting)
%getParamIndices Get the Indices for All Moduli and Characteristic Times
%   This function creates a list of the indices for the Moduli and
%   Characteristic times inside an array of size beta with the settings
%   dictated by elasticSetting and fluidSetting (which are 'y' and 'n'
%   depending on whether an elastic term and steady-state fluidity are
%   included in the viscoelastic model).

if strcmp(elasticSetting,'y')
    if strcmp(fluidSetting,'y')
        tauInds = (3:2:length(beta)-1);
        modulusInds = horzcat(1,(2:2:length(beta)-1));
    else
        tauInds = (3:2:length(beta));
        modulusInds = horzcat(1,(2:2:length(beta)));
    end
else
    if strcmp(fluidSetting,'y')
        tauInds = (2:2:length(beta)-1);
        modulusInds = (1:2:length(beta)-1);
    else
        tauInds = (2:2:length(beta));
        modulusInds = (1:2:length(beta));
    end
end

end

