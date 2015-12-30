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

R   = 100;   % number of runs

dt = T/N;

fprintf(' Longstaff-Schwartz Monte Carlo method \n');
P = zeros(R, 1);
for i = 1:1:R
    fprintf('\nRun No. %d:\n', i);
    P(i) = single_run_in_the_money(N, M, dt, S0, r, sig, K);
end

val = mean(P);
std_err = std(P) / sqrt(length(P)); % standard error of the mean

fprintf('Average of %d runs: %f\n', R, val);
fprintf('Confidence interval (2 standard errors of mean): [%f, %f]\n', ...
        val - 2 * std_err, val + 2 * std_err);
