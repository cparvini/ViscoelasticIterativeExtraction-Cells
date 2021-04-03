function [x_log] = log_scale(x,t,tr,st)
%log_scale(variable, time array, dt, stopping time)
%   This function takes an input of two arrays, the first being the
%   variable of interest and the second being the time array for that
%   variable. Using a while loop, it resamples the arrays into logarithmic
%   form. The final setting chooses whether or not to skip negative values
    
    n_log = length(x);
    n_logMax = 10;
    i = 1;
    j = 1;
    x_log = [];
    
    while i < n_log
        % How many points do we have in this decade?
        magPoints = (find( floor(log10(t)) == floor(log10(j*tr)) ));
        if length(magPoints) < 10
            stepsize = 1;
        else
            stepsize = median( floor( diff( linspace(magPoints(1), magPoints(end), 10) ) ) );
        end
        
        % If you are between the first and last time...
        if t(i) >= j*tr && t(i) <= st
            
            x_log = [x_log x(i)];
            j = j + 1;
            
            % If you've hit the top of this log range...
            if j == n_logMax
                tr = tr * 10;
                j = 1;
            end
        end

        i = i + stepsize;
        
    end

end

