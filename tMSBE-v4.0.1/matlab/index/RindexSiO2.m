function [n] = RindexSiO2(x)
n = sqrt( 1.28604141 + 1.07044083.*power(x,2)./(power(x,2)-1.00585997e-2) + 1.10202242.*power(x,2)./(power(x,2)-100));
end

