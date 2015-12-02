
%
% modelling constants
%

close all;
clear all;

r    = 0.05;
sig  = 0.2;
T    = 1;
X    = 1;
Smax = 2;

%
% grid
%

J = 201;
N = 4001;

% J = 401;     % uncomment to run with twice as many grid points and
% N = 16001;   % four times as many timesteps, for numerical stability

S  = linspace(0,Smax,J)';
dS = S(2)-S(1);
t  = linspace(0,T,N);
dt = t(2)-t(1);

j  = 1:J-1;
jm = max(j-1,1);
jp = j+1;

V0 = max(X-S,0);
V  = V0*ones(1,N);
Sb = zeros(1,N);

%
% time-marching
%

Sb(N) = 1;

for n = N:-1:2
  V(j,n-1) = (1-dt*r)*V(j,n)                ...
           + 0.5*dt*r*S(j).*(V(jp,n)-V(jm,n))/dS ...
           + 0.5*dt*sig^2*S(j).^2.*(V(jp,n)-2*V(j,n)+V(jm,n))/dS^2; 
  V(j,n-1) = max(V0(j),V(j,n-1));

  Sb(n-1) = max(find(V(:,n-1) < X-S+1e-10)) * dS;
end

%
% plot results
%

N1 = 1;
N2 = 1+(N-1)/4;
N3 = 1+(N-1)/2;
N4 = 1+3*(N-1)/4;
N5 = N;

figure
plot(S,V(:,N1), S,V(:,N2), S,V(:,N3), S,V(:,N4), S,V(:,N5))
axis([0.8 1.2 0 0.2])
title('American put option')

%figure
%plot(t,Sb,t,Sb-dS)
%xlabel('t'); ylabel('S_b')
%axis([0 1 0.8 1])
%title('American put option')

fprintf(' explicit finite difference method \n');
fprintf(' option value at S0=1 is %f \n',V(1+(J-1)/2,1));
