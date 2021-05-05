function [n] = RindexBaTiO3(x)
n = sqrt( 1 + 4.187.*power(x,2)./(power(x,2)-power(0.223,2)));
end

