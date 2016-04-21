function [r,g,b] = jetplot(Coeff)

if Coeff <= 8
    r = 0;
    g = 0;
    b = 0.5 + 0.5/8*Coeff;
elseif ( Coeff > 8 ) && ( Coeff <= 24 )
    r = 0;
    g = 1/(24-8)*(Coeff-8);
    b = 1;
elseif ( Coeff > 24 ) && ( Coeff <= 40 )
    r = 1/(40-24)*(Coeff-24);
    g = 1;
    b = 1 - r;
elseif ( Coeff > 40 ) && ( Coeff <= 56 )
    r = 1;
    g = 1 - 1/(56-40)*(Coeff-40);
    b = 0;
else
    r = 1 - 1/(64-56)*(Coeff-56);
    g = 0;
    b = 0;
end

