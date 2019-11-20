function [q] = VVC_corrected(V,qmin,qmax,Vmin,Vmax)

if V < Vmin
    q = -qmax;
elseif Vmin <= V && V <= Vmax
    q = (qmax - qmin)/(Vmax - Vmin)*(V - Vmin) + qmin;
elseif V > Vmax
    q = qmax;
end

