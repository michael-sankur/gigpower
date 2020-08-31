function [uk, rhok, sigmak, xik, uhatk, ak, Jdotk, uhatdotk] = ...
    esfunction(t, dt, hpf, lpf, f, ...
    kint, akm1, Jk, Jkm1, Jdotkm1, ...
    rhokm1, sigmakm1, xikm1, uhatkm1, a0, ...
    theta, sw, beta)

% calculate angular frequency:
w = 2*pi*f;

% calculate derivative of objective function
Jdotk = (2/dt)*(Jk - Jkm1) - Jdotkm1;

% extract the effect of the probing signal in the objective function
% do this by passing the signal through a highpass filter
% rhok = (Jk - Jkm1 - (hpf*dt/2-1)*rhokm1)/(1+hpf*dt/2)
rhok = (2*(Jk - Jkm1) + (2 - hpf*dt)*rhokm1)/(2 + hpf*dt);

% the resulting signal is a sinusoid, multiply by a sinusoid of the same frequency
% this results in a cos**2 term, that has a DC component (we call this demodulation)
sigmak = rhok*cos(w*t + theta);

% pass the demodulated signal through a lowpass filter, to eliminate noise and "jumpiness"
xik = (dt*lpf*(sigmak + sigmakm1) - (dt*lpf - 2)*xikm1)/(2 + dt*lpf);

% pass the resulting signal through an integrator - this approximates a gradient descent
uhatk = uhatkm1 + kint*dt/2*(xik + xikm1);

uhatdotk = kint*xik;

if sw == 0
    ak = a0;
elseif sw == 1
    ak = akm1*(1 - beta*dt);
    ak = akm1*(2 - beta*dt)/(2 + beta*dt);
end

% modulation - add the perturbation to the next control setpoint
uk = (uhatk + ak*cos(w*t + theta));


end
