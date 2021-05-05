function [n] = RindexGaAs(lambda,T)
T0 = 300;
n = sqrt(3.5 + 7.4969*power(lambda,2)./(power(lambda,2)-power(0.4082,2)) + 1.9347*power(lambda,2)./(power(lambda,2)-power(37.17,2)) ) + 2.67e-4.*(T-T0) ;

end

