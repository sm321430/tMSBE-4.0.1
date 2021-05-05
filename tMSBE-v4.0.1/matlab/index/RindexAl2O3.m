function [n] = RindexAl2O3(x)
n = sqrt(1+1.4313493*power(x,2)./(power(x,2)-power(0.0726631,2))+0.65054713*power(x,2)./(power(x,2)-power(0.1193242,2))+5.3414021*power(x,2)./(power(x,2)-power(18.028251,2)));
end

