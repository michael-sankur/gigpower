function [Xavg] = avgfunction(kt, t, dt, f, X)

w = 2*pi*f;
tau = 2*pi/w;

if kt == 1
    Xavg = dt*X(1);
elseif kt == 2
    Xavg = dt/2*(X(1) + X(2));
elseif kt >= 3 && kt*dt <= tau
    Xavg = dt*(1/2)*(X(1) + 2*sum(X(2:kt-1)) + X(kt));
elseif kt >= 3 && kt*dt > tau
%     f
%     kt
%     kt-floor(tau/dt)
    Xavg = f*dt*(1/2)*(X(kt-floor(tau/dt)) + 2*sum(X(kt-floor(tau/dt)+1:kt)) + X(kt));
end

end