clear; clc;


k_sqrd = 156*(2^(-7.0/3.0)) - 42*2^( -4.0/3.0);


% Coarse Grid parameters
a = 0;
b = 10000;
dt = 0.2;
t = a:dt:b;
k = sqrt(k_sqrd);

% Memory Function Parameters
epsilon = 10^(-16); % Error Tolerance
r = 0.99; % Decay Rate
n = 0;
T = 2^n;
Tc = T * sqrt(log(r)/log(epsilon)); % r-life decay time

% Initialize Memory function on coarse Grid
x = theta(t,k);
%x = trunc_theta(t,k,r,Tc);
x(1) = 1;

%plot(t,x)
N = length(x);
L = floor(N/2);
M = 100;

Y = zeros(N-L,L+1);

for jj = 1:N-L
    Y(jj,:) = x(jj:(jj+L));
    
end

[~, S, ~] = svds(Y,M);
Y = NaN;


save("soe_approx_theta")

function out = theta(t,k)
 out = besselj(1,2*k.*t)./(k.*t);
end



