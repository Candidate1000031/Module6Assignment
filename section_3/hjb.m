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
M = 10;
% number of mesh intervals
N = 9;
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
    a = (sigma * p(j) * W) .^ 2 / 2;
    b = pi + W * (r + p(j) * sigma * xi);

    alpha_central = a / (dW ^ 2) - b / (2 * dW);
    beta_central = a / (dW ^ 2) + b / (2 * dW);

    alpha_fb = a / (dW ^ 2) + max(0, - b / (2 * dW));
    beta_fb = a / (dW ^ 2) + max(0, b / (2 * dW));

    % find indices where central differences are both positive
    mask = (alpha_central >= 0) & (beta_central >= 0);
    inverse_mask = ones(size(W)) - mask;

    alpha = alpha_central .* mask + alpha_fb .* inverse_mask;
    beta = beta_central .* mask + beta_fb .* inverse_mask;

    % TODO possible take time into account already here
    gamma_ = -alpha - beta;

    coeff_matrices{j} = ...
        spdiags([alpha gamma_ beta], -1:1, N + 1, N + 1);
end

% incorporate the boundary conditions
