clear; clc;


k_sqrd = 156*(2^(-7.0/3.0)) - 42*2^( -4.0/3.0);


% Coarse Grid parameters
a = 0;
b = 10;
dt = 0.02;
t = a:dt:b;
k = sqrt(k_sqrd);

% Memory Function Parameters
epsilon = 10^(-16); % Error Tolerance
r = 0.99; % Decay Rate
T = 2^15;
Tc = T * sqrt(log(r)/log(epsilon)); % r-life decay time

% Initialize Memory function on coarse Grid
%x = theta(t,k);
x = trunc_theta(t,k,r,Tc);
x(1) = 1;

plot(t,x)


function out = theta(t,k)
 out = besselj(1,2*k.*t)./(k.*t);
end

function out = trunc_theta(t,k,r,Tc)
    out = r.^((t./Tc).^2).*theta(t,k);
end