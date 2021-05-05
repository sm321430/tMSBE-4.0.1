function [n] = RindexDiam(x)

n = sqrt( 1 + 4.3356.*power(x,2)./(power(x,2)-power(0.1060,2)) + 0.3306.*power(x,2)./(power(x,2)-power(0.1750,2)) );
end

