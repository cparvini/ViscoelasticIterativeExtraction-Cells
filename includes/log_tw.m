function [omega] = log_tw(de0, maxi)
%log_tw( time vector, time resolution )
%   This function generates a frequency or time vector in log scale.

    omega = [];
    w = [];
    de = [];
    prnts = [];
    
    nn = 20;
    w = de0;
    de = de0;
    prints = 1;

    omega = [omega;de0];

    while w < maxi
        w = w + de;
        if w < maxi
            omega = [omega;w];          
            prints = prints + 1;
        end
        if prints == nn
            de = de * 10;
            prints = 1;
        end
    end

end

