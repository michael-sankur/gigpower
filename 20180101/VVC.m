function [q] = VVC(V,qmax,Vmin,Vmax)

if V <= Vmin
    q = -qmax;
elseif V > Vmin <= Vmax
    q = (V-1)*(2*qmax)/(Vmax - Vmin);
elseif V > Vmax
    q = qmax;
end

