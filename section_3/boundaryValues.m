function [V_Wmax] = boundaryValues(pi, r, gamma, Wmax, t)
%BOUNDARYVALUES calculates the boundary values for W = Wmax
    tau = T - t;

    c = 2 * pi / r;
    alpha = exp(2 * r * tau);
    beta = -(gamma + c) * exp(r * tau) + c * exp(2 * r * tau);
    delta = -((pi * (gamma + c)) / r) * (exp(r * tau) - 1) + ...
        ((pi * c) / (2 * r)) * (exp(2 * r * tau) - 1) + gamma ^ 2 / 4;

    V_Wmax = alpha * Wmax ^ 2 + beta * Wmax + delta;
end
