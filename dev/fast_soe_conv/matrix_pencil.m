clear; clc;


k_sqrd = 156*(2^(-7.0/3.0)) - 42*2^( -4.0/3.0);


% Coarse Grid parameters
a = 0;
b = 10000;
dt = 0.5;
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

%Y2 = zeros(N-L,L);
%Y1 = zeros(N-L,L);
Y = zeros(N-L,L+1);

for jj = 1:N-L
    Y(jj,:) = x(jj:(jj+L));
    %Y1(jj,:) = x(jj:(jj + L-1));
    %Y2(jj,:) = x((jj+1):(jj + L));
    
end

[U, S, V] = svds(Y,M);
%S(M,M)
%%
U = NaN;
Y = NaN;

V1 = V(1:(end-1),:);
V2 = V(2:end,:);

V = NaN;


%E = eig(pinv(V1,S(M,M))*V2);
E = eig(pinv(V1,S(M,M))*V2);
V1 = NaN;
V2 = NaN;
%S = NaN;
E = E';
%Y1 = U*S*V1';
%Y2 = U*S*V2';
%E = eig(pinv(Y1)*(Y2));
%E = E(1:M);

Z = ones(N,M);
for jj = 2:N
   Z(jj,:) = E.^(jj-1);
end

E = NaN;

R = Z\x';

Y_out = zeros(N,1);
for jj = 1:N
   for ii = 1:M
       Y_out(jj) = Y_out(jj) + R(ii).*Z(2,ii).^(jj-1);
   end
end

Z_vec = Z(2,:);
Z = NaN;


Y_out = real(Y_out);

%save("soe_approx_theta_T_is" + num2str(T))
save("soe_approx_theta")
%%
figure
subplot(2,3,1)
plot(t, (Y_out - x')./x')
title('Relative error')


subplot(2,3,2)
plot(t,abs(Y_out - x'))
title('Absolute Error')

subplot(2,3,3)
plot(t,Y_out)
title('Sum of Exp Approximation')

subplot(2,3,4)
plot(t,x')
title('True Samples')


% Initialize fine grid
a2 = a;
b2 = b;
dt2 = 0.001;
t2 = a2:dt2:b2;

% Initialize memory on fine grid
x2 = theta(t2,k);
%x2 = trunc_theta(t2,k,r,Tc);
x2(1) = 1;




subplot(2,3,5)
SoE = SoE_approx(t2,R,Z_vec,dt);
plot(t2,x2 - SoE_approx(t2,R,Z_vec,dt))
title('Cont SOE approx Error')

function out = SoE_approx(t,R,Z,new_dt)
 out = zeros(size(t));
 s = log(Z);
 parfor jj = 1:length(out)
       out(jj) = dot(R',exp( t(jj) .*s./new_dt));
 end
out = real(out);
end

function out = theta(t,k)
 out = besselj(1,2*k.*t)./(k.*t);
end

function out = trunc_theta(t,k,r,Tc)
    out = r.^((t./Tc).^2).*theta(t,k);
end


