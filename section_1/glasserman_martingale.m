function M_k = glasserman_martingale(M_k_minus_1, X_k_minus_1, X_k, beta_k, ...
                                     dt, r, sig, K, P)
% GLASSERMAN_MARTINGALE computes the values of the Glasserman martingale at the
% k-th time point.
%
% M_k_minus_1 - The values of the Glasserman martingale at the k-1-th time point.
% X_k_minus_1 - The values of the underlying at time k-1. 
% X_k - The values of the underlying at time k. 
% beta_k - The time k regression coefficients from Longstaff-Schwartz.
% dt - Length of time step.
% r - Risk-free rate.
% sig - Volatility of underlying asset.
% K - Strike.
% P - Number of mini-paths.

    get_continuation_value = ...
        @(X) [ones(size(X)), 1 - X,  1 / 2 * (2 - 4 * X - X .^ 2)] * beta_k;

    get_V = @(X) max(K - X, get_continuation_value(X));

    % values from mini-paths
    mini_paths = get_mini_paths(X_k_minus_1, dt, r, sig, P);

    for i = P:-1:1
        X = mini_paths(:, i);
        V_k_p(:, i) = get_V(X);
    end

    M_k = M_k_minus_1 + get_V(X_k) - sum(V_k_p, 2) / P;
end

function X_k = get_mini_paths(X_k_minus_1, dt, r, sig, P)
% GET_MINI_PATHS - Get P mini paths starting from X.
%
% X_k_minus_1 - The values of the underlying at the starting point. 
% dt - Length of time step.
% r - Risk-free rate.
% sig - Volatility of underlying asset.
% P - Number of mini-paths.

    M = length(X_k_minus_1);

    % Brownian increments
    dW = sqrt(dt) * randn(M, P);                              

    % mini-paths
    X_k = repmat(X_k_minus_1, 1, P) .* exp((r - sig ^ 2 / 2) * dt + sig * dW);
end
