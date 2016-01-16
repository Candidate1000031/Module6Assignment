function [Vnew, k] = policyIteration(coeff_matrices, rhs, Vini)
    n = length(Vini);

    Vnew = Vini;
    Vold = 2 * Vini;
    
    eps = 1e-6;
    maxit = 50;
    it = 0;

    J = length(coeff_matrices);

    while (norm(Vold - Vnew) / norm(Vini) > eps && it < maxit)
        Vold = Vnew;
        
        % k will hold the index of the p_k that minimizes the lhs
        k = ones(n, 1);

        % Vmin will hold the minimal value of the lhs
        Vmin = coeff_matrices{1} * Vold;

        % this will hold the coefficients for those p_k that minimize the rhs
        mat = coeff_matrices{1};

        for j = 2:J
            V = coeff_matrices{j} * Vold;

            mask = V < Vmin;
            inverse_mask = V >= Vmin;
            matmask = diag(mask);
            inverse_matmask = diag(inverse_mask);

            k = k .* inverse_mask + j * ones(n, 1) .* mask;
            Vmin = Vmin .* inverse_mask +  V .* mask;

            mat = inverse_matmask * mat + matmask * coeff_matrices{j};
        end

        Vnew = mat \ rhs;
        it = it + 1;
    end
    fprintf('Iterations needed: %d\n', it);
end
