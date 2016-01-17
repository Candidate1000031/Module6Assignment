% Solve HJB equation, using central differences "as much as possible"
% Parameters:
r = .03;
sigma = .15;
xi = .33;
pi = .1;
W0 = 1;
T = 20;
gamma = 14.47;
pmax = 1.5;
Wmax = 5;

% number of timesteps
M = 1600;
% number of mesh intervals
N = 100;
% number of discrete values for p
J = 8;

% length of timestep
dt = T / M;
% length of mesh interval
dW = Wmax / N;

% discretize W, t, p
W = linspace(0, Wmax, N + 1)';
t = linspace(0, T, M + 1);
p = linspace(0, pmax, J);

% initialize V
V = initializeV(gamma, W);

% cells for storing the coefficient matrices
coeff_matrices = cell(1, J);

% parameters and coefficients
for j = 1:J
    a = ((sigma * p(j) * W) .^ 2) / 2;
    b = pi + W * (r + p(j) * sigma * xi);

    alpha_central = a / (dW ^ 2) + b / (2 * dW);
    beta_central = a / (dW ^ 2) - b / (2 * dW);

    alpha_fb = a / (dW ^ 2) + max(0, b / dW);
    beta_fb = a / (dW ^ 2) + max(0, -b / dW);

    % find indices where central differences are both positive
    mask = (alpha_central >= 0) & (beta_central >= 0);
    inverse_mask = ones(size(W)) - mask;

    alpha = alpha_central .* mask + alpha_fb .* inverse_mask;
    beta = beta_central .* mask + beta_fb .* inverse_mask;
    gamma_ = -alpha - beta - (1 / dt);

    coeff_matrices{j} = ...
        spdiags([alpha gamma_ beta], -1:1, N + 1, N + 1)';

    % last row only for enforcing boundary condition
    coeff_matrices{j}(N + 1, :) = 0;
    coeff_matrices{j}(N + 1, N + 1) = 1;
end


V_Wmax = boundaryValues(pi, r, gamma, Wmax, t, T);

% loop over time steps
for m = 1:M
    rhs = -V / dt;
    rhs(N + 1) = V_Wmax(M + 1 - m);

    [V, k] = policyIteration(coeff_matrices, rhs, V);
end

figure
plot(W, V);
xlabel('Wealth W');
ylabel('V(W, 0)');
figure
plot(W, k * pmax / (J - 1))
xlabel('Wealth W');
ylabel('Optimal control p(W, t = 0)');
