function [uhatk,vhatk] = rectify_what_generation(tk,uk,vk,uhatk,vhatk,fes,a0)

if uk >= 0 && 2*pi*fes*tk ~= pi/2 && 2*pi*fes*tk ~= -pi/2
    hk = vhatk - uhatk*tan(2*pi*fes*tk);
%     ck = a0*(cos(2*pi*fes*tk) + 1j*sin(2*pi*fes*tk));
    ck = a0*exp(1j*2*pi*fes*tk);
    whatk = 1j*hk - ck;
    uhatk = real(whatk);
    vhatk = imag(whatk);
end

end

