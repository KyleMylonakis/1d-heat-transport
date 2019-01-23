clear; clc;


k_sqrd = 156*(2^(-7.0/3.0)) - 42*2^( -4.0/3.0);


% Coarse Grid parameters
a = 0;
b = 10000;
dt = 0.05;
t = a:dt:b;
k = sqrt(k_sqrd);

% Memory Function Parameters

epsilon = 10^(-16); % Error Tolerance
r = 0.99; % Decay Rate
n = 0;
T = 2^n;
Tc = T * sqrt(log(r)/log(epsilon)); % r-life decay time

% Initialize Memory function on coarse Grid
%x = theta(t,k);
x = trunc_theta(t,k,r,Tc);
x(1) = 1;

%plot(t,x)
N = length(x);
L = floor(N/2);
M = 55;

%Y2 = zeros(N-L,L);
%Y1 = zeros(N-L,L);
Y = zeros(N-L,L+1);

for jj = 1:N-L
    Y(jj,:) = x(jj:(jj+L));
    %Y1(jj,:) = x(jj:(jj + L-1));
    %Y2(jj,:) = x((jj+1):(jj + L));

end

[~, S, V] = svds(Y,M);
%S(M,M)
%%
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

gamma = log(Z_vec)/dt;

gamma_real = real(gamma);
gamma_imag = imag(gamma);
R_real = real(R);
R_imag = imag(R);

gamma_real_fileID = fopen('gamma_real.bin', 'w');
fwrite(gamma_real_fileID,gamma_real,'double');
fclose(gamma_real_fileID);

gamma_imag_fileID = fopen('gamma_imag.bin', 'w');
fwrite(gamma_imag_fileID,gamma_imag,'double');
fclose(gamma_imag_fileID);

R_real_fileID = fopen('R_real.bin', 'w');
fwrite(R_real_fileID,R_real,'double');
fclose(R_real_fileID);

R_imag_fileID = fopen('R_imag.bin', 'w');
fwrite(R_imag_fileID,R_imag,'double');
fclose(R_imag_fileID);

save("soe_approx_theta_T_is_" + num2str(T))
%save("soe_approx_theta")

function out = theta(t,k)
 out = besselj(1,2*k.*t)./(k.*t);
end

function out = trunc_theta(t,k,r,Tc)
    out = r.^((t./Tc).^2).*theta(t,k);
end


