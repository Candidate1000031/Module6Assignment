function [V] = initializeV(gamma, W)
%INITIALIZEV initializes values for V at time T
    V = (W - gamma / 2) .^ 2;
end
