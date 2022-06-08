function [theta, U] = possibilistic_2(X, m, theta_init, eta, max_iter)
%   Summary of this function goes here
%   X: LxN matrix
%   theta_init: mxL matrix
%   eta: 1xm matrix
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
    U = exp(- (d ./ (ones(N,1)*eta)));
    
    
    % compute thetas
    for j=1:m
        A = sum( (U(:,j)*ones(1,L))' .* X, 2);
        theta(j,:) = A' ./ sum(U(:,j));
    end
    
    e = norm(theta_old - theta, inf);
    iter = iter+1;
    
    iter, e
end