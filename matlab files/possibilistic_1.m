function [theta, U] = possibilistic_1(X, m, theta_init, eta, q, max_iter)
%   Summary of this function goes here
%   X: LxN matrix
%   theta_init: mxL matrix
%   eta: 1xm matrix
%   q>1
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
    U = (ones(N,m) + (d ./ (ones(N,1)*eta)) .^ (1 /(q-1)) ) .^(-1);
    
    % compute thetas
    U_q = U .^ q;
    for j=1:m
        A = sum( (U_q(:,j)*ones(1,L))' .* X, 2);
        theta(j,:) = A' ./ sum(U_q(:,j));
    end
    
    e = norm(theta_old - theta, inf);
    iter = iter+1;
    
    iter, e
end