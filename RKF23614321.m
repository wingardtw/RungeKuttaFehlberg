function [ Tout, Xout, DXout, info ] = RKF23614321( T0, Tfinal, X0, DX0, tol, A, Mu, omega )
%RKF23614321 Solves van der Pol equation using Runge-Kutta Felberg
%technique
%   T0 - initial time
%   Tfinal - final time
%   X0 - value of x at time T0
%   DX0 - value of x' at time T0
%   tol - error tolerance
%   A,Mu,omega - Variables in van der Pol equation
info = 'nothanks';
assert(Mu >= 0, 'Mu must be positive or zero')
t = T0;
u = [X0, DX0];
Tout = t;
Xout = u(1);
DXout = u(2);
hmax = (Tfinal - T0)/1000; %Maximum stepsize corresponds to 1000 steps
hmin = 16*eps(t) + eps;
h = hmax;
FLAG = 1;

% Full equation -- Mu > 0, A != 0
while FLAG == 1
    K1 = h*vdp(t,u,Mu,A,omega);
    K2 = h*vdp(t+(1/4)*h, u + (1/4)*K1,Mu,A,omega);
    K3 = h*vdp(t + (3/8)*h, u + (3/32)*K1 + (9/32)*K2, Mu, A, omega);
    K4 = h*vdp(t + (12/13)*h, u + (1932/2197)*K1 - (7200/2197)*K2 + (7296/2197)*K3, Mu, A, omega);
    K5 = h*vdp(t + h, u + (439/216)*K1 - 8*K2 + (3680/513)*K3 - (845/4104)*K4, Mu, A, omega);
    K6 = h*vdp(t + (1/2)*h, u - (8/27)*K1 + 2*K2 - (3544/2565)*K3 + (1859/4104)*K4 - (11/40)*K5, Mu, A, omega);
    R = (1/h)*abs(K1/360 - (128/4275)*K3 - (2197/75240)*K4 +(1/50)*K5 + (2/55)*K6);
    
    %Accepted
    if R <= tol
        t = t + h;
        u = (u + (25/216)*K1 + (1408/2565)*K3 + (2197/4104)*K4 - K5/5);
        Tout = [Tout t];
        Xout = [Xout u(1)];
        DXout = [DXout u(2)];
    end
    
    %New h
    del = 0.84*nthroot([tol,tol]/R, 4);
    if del <= 0.1
        h = 0.1*h;
    elseif del >= 4
        h = 4*h;
    else
        h = del*h;
    end
    
    if h > hmax
        h = hmax;
    end
    
    if t >= Tfinal
        FLAG = 0;
    elseif (t+h > Tfinal)
        h = Tfinal - t;
    elseif h < hmin
        FLAG = 0;
        info = 'something wrong';
        return
    end
    
end
end


%Auxiliary Functions

function [du] = vdp(t,u,Mu,A,omega)
du = [u(2),Mu*(1-u(1)^2)*u(2) - u(1) + A*sin(omega*t)];
end