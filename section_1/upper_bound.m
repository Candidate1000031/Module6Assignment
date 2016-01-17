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
P   = 100;   % number of mini paths

dt = T/N;

%
% Phase 1: run Longstaff-Schwartz to get regression parameters for approximate
% continuation values
%
fprintf('Running Longstaff-Schwartz lower bound calculations\n');
[~, betas] = single_run(N, M, dt, S0, r, sig, K);


%
% Phase 2: calculate Glassermann martingale approximation
%
% paths
dW = sqrt(dt)*randn(N,M);                              % Brownian increments
S  = cumprod([repmat(S0,1,M); exp((r-sig^2/2)*dt+sig*dW)]);  % paths

mart = zeros(N, M);
for k = 2:1:N
    mart(k, :) = glasserman_martingale(mart(k - 1, :)', S(k - 1, :)', ...
                                       S(k, :)', betas(k, :)', dt, r, sig, ...
                                       K, P);
end

%
% Phase 3: estimate expectation using mean
%
res = max(max(K - S(1:N, :), 0) - mart);
val = mean(res);
sd = std(res);


fprintf('Upper bound estimate according to Glasserman \n');
fprintf('option value = %f, std dev = %f \n', val, sd);
