%
% Based on code originally downloaded from 
% http://www.mathworks.com/matlabcentral/fileexchange/
%        16476-pricing-american-options/content/AmericanOptLSM.m
%
% AmericanPutLSM - Price an american put option via Longstaff-Schwartz Method
%
% Inputs:
%
%   S0      Initial asset price
%   K       Strike Price
%   r       Interest rate
%   T       Time to maturity of option
%   sigma   Volatility of underlying asset
%   N       Number of timesteps
%   M       Number of paths 

S0  = 1.0;   % initial asset price
K   = 1.0;   % strike
r   = 0.05;  % risk-free interest rate
T   = 1.0;   % time to maturity
sig = 0.2;   % volatility

N   = 64;    % number of timesteps
M   = 1e5;   % number of paths

dt = T/N;

%
% Phase 1: least squares regression
%

% generate paths

dW = sqrt(dt)*randn(N,M);                              % Brownian increments
S  = cumprod([repmat(S0,1,M); exp((r-sig^2/2)*dt+sig*dW)]);  % paths

P = max(K-S(end,:),0);  % put option payoff

% loop backwards in time

for n = N:-1:2
  P    = exp(-r*dt)*P;                            % discount 
  ind  = 1:M;                                     % use all paths
  X    = S(n,ind)';
  Y    = P(ind)';
  A    = [ ones(size(X)) (1-X) 1/2*(2-4*X-X.^2)]; % 3 basis functions
  beta = A\Y;                                     % linear regression
  betas(n,:) = beta';
  C    = A*beta;                                  % continuation value

  ind  = ind(find(K-X > C));
  P(ind) = K-S(n,ind);
end

P   = exp(-r*dt)*P;
val = sum(P)/M;
sd  = sqrt((sum(P.^2)/M - val^2)/M);

fprintf(' Longstaff-Schwartz Monte Carlo method \n');
fprintf(' phase 1: option value = %f, std dev = %f \n',val, sd);

