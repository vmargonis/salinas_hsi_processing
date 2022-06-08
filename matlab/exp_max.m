function [theta, U] = exp_max(X, m, theta_init, sigma, P, max_iter)
% X: LxN data
% theta: mxL, initial means
% sigma: mx1, initial co-variances
% P: mx1, initial prior probabilities
[L,N]=size(X);


% preallocate U, d, B
U = zeros(N,m);
B = zeros(N,m);
d = zeros(N,m);

% create 3-dim matrix for covariance arrays
% S will be an LxLxm array
S = sigma(1,1)^2 .* eye(L);
for j=2:m
    S(:,:,j) = sigma(j,1) .* eye(L);
end

e=1;
iter=0;
theta = theta_init;

while(e > 0.001) && (iter < max_iter)
    theta_old = theta;
    P_old = P;
    
    % compute mahalanobis distances from representatives   
    for j=1:m
        y = X - theta(j,:)'*ones(1,N);
        a = S(:,:,j) \ y;
        for i=1:N
            d(i,j) = y(:,i)' * a(:,i);
        end
    end    
    
    d = d + (d==0)*10^(-10);
    
    % compute posterior probabilities (expectation step)
    % -----------------------------------------------------------
    for j=1:m
        B(:,j) = det(S(:,:,j))^(-0.5) * P(j,1) * exp(-0.5 * d(:,j));
    end
    
    U = B ./ (sum(B,2) * ones(1,m));
    % -----------------------------------------------------------
    
    % maximization step
    % -----------------------------------------------------------
    gamma = sum(U);
    for j=1:m
        % compute new means
        m1 = (U(:,j) * ones(1,L)) .* X';
        theta(j,:) = sum(m1) ./ gamma(1,j);
        
        % compute new covariances
        m2 = zeros(L,L);
        for i=1:N
            xi_mj = X(:,i) - theta(j,:)'; % Lx1
            m2 = m2 + U(i,j) .* (xi_mj * xi_mj'); %LxL
        end
        S(:,:,j) = m2 ./ gamma(1,j);
    end
    
    % compute new prior probabilities
    P = 1/N * gamma';
    % -----------------------------------------------------------

    e = norm(theta_old - theta) + norm(P_old - P);
    iter = iter+1;
    
    e, iter
end

end
