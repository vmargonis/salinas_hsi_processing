function [theta, U] = fuzzy(X, m, theta_init, q, max_iter)
%   Summary of this function goes here
%   X: LxN matrix
%   theta_init: mxL matrix

[L,N]=size(X);

% preallocate U, d
U = zeros(N,m);
d = zeros(N,m);

e=1;
iter=0;
theta = theta_init;

while(e > 0.001) && (iter < max_iter)
    theta_old = theta;
    
    % compute distances from representatives
    for j=1:m
        d(:,j) = sum( (X - theta(j,:)'*ones(1,N)).^2 )' ;
    end
    d = d + (d==0)*10^(-10);
    
    % compute U
    for j=1:m
        D = ((d(:,j)*ones(1,m)) .* (d .^ (-1))) .^ (1/(q-1));
        U(:,j) = ones(N,1) ./ sum(D,2);
    end
    
    % update theta
    U_q = U .^ q;
    for j=1:m
        A = sum( (U_q(:,j)*ones(1,L))' .* X, 2);
        theta(j,:) = A' ./ sum(U_q(:,j));
    end
    
    e = norm(theta_old - theta, inf);
    iter = iter+1;
    iter, e
end

end