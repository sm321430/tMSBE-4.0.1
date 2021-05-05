function [n] = RindexSi3N4(x)
n = sqrt(1 + 2.8939*power(x,2)./(power(x,2)-power(139.67e-3,2)));
end

