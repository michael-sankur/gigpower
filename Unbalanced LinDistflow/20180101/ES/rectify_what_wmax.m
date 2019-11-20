function [uhatk,vhatk] = rectify_what_wmax(tk,uk,vk,uhatk,vhatk,fes,a0,wmaxpu)

if sqrt(uk^2 + vk^2) >= wmaxpu
    wk = uk + 1j*vk;
    %     whatk = uhatk + 1j*vhatk;
    %     ck = wk - whatk;
    %     ck = a0*(cos(2*pi*fes(ph,kn)*time(kt)) + 1j*sin(2*pi*fes(ph,kn)*time(kt)));
    ck = a0*exp(1j*2*pi*fes*tk);
    whatk = wmaxpu*wk/abs(wk) - ck;
    uhatk = real(whatk);
    vhatk = imag(whatk);
end

end

