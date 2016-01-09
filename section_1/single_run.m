function [val, betas]  = single_run(N, M, dt, S0, r, sig, K)
% SINGLE_RUN computes Longstaff-Schwartz estimate of an American call option
% for a single run.

    %
    % Phase 1: least squares regression
    %
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

        ind  = ind(K-X > C);
        P(ind) = K-S(n,ind);
    end

    P   = exp(-r*dt)*P;
    val = sum(P)/M;
    sd  = sqrt((sum(P.^2)/M - val^2)/M);

    fprintf(' Phase 1: option value = %f, std dev = %f \n',val, sd);

    %
    % Phase 2: final price evaluation
    %
    % second set of paths for price evaluation
    dW = sqrt(dt)*randn(N,M);                              % Brownian increments
    S  = cumprod([repmat(S0,1,M); exp((r-sig^2/2)*dt+sig*dW)]);  % paths

    P = max(K-S(end,:),0);  % put option payoff

    % loop backwards in time
    for n = N:-1:2
        P    = exp(-r*dt)*P;                            % discount 
        ind  = 1:M;                                     % use all paths
        X    = S(n,ind)';
        A    = [ ones(size(X)) (1-X) 1/2*(2-4*X-X.^2)]; % 3 basis functions
        beta = betas(n,:)';
        C    = A*beta;                                  % continuation value

        ind  = ind(K-X > C);
        P(ind) = K-S(n,ind);
    end

    P   = exp(-r*dt)*P;
    val = sum(P)/M;
    sd  = sqrt((sum(P.^2)/M - val^2)/M);

    fprintf(' Phase 2: option value = %f, std dev = %f \n',val, sd);
end
