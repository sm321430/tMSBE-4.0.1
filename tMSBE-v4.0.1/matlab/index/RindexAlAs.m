function [n] = RindexAlAs(x,T)
T0 = 300;
n = sqrt( 2.0792 + 6.0840*power(x,2)./(power(x,2)-power(0.2822,2)) + 1.900*power(x,2)./(power(x,2)-power(27.62,2)) ) + 1.43e-4.*(T-T0);
end

