function [dydx, d2ydx2] = numerical_derivatives(Y, x)
%NUMERICAL_DERIVATIVES Computes numerical derivatives of a matrix Y w.r.t a vector x
%   Y - Matrix with column data (n x m), n is number of points, and m is number of variables
%   x - vector for differentiation (n x 1)

% Ensure x is a column vector
x = x(:);
n = length(x);
m = size(Y, 2); % Number of columns in Y

% init
dydx = zeros(n, m);
d2ydx2 = zeros(n, m);

for i = 1:m
    y = Y(:, i);
    % Central difference for first and second derivatives
    dydx(2:n-1, i) = (y(3:n) - y(1:n-2)) ./ (x(3:n) - x(1:n-2));
    d2ydx2(2:n-1, i) = (y(3:n) - 2*y(2:n-1) + y(1:n-2)) ./ ((x(3:n) - x(1:n-2)) / 2).^2;
    
    % Forward difference for first point to keep dimentions
    dydx(1, i) = (y(2) - y(1)) / (x(2) - x(1));
    d2ydx2(1, i) = (y(3) - 2*y(2) + y(1)) / ((x(3) - x(1)) / 2)^2;
    
    % Backward difference for last point to keep dimentions
    dydx(n, i) = (y(n) - y(n-1)) / (x(n) - x(n-1));
    d2ydx2(n, i) = (y(n) - 2*y(n-1) + y(n-2)) / ((x(n) - x(n-2)) / 2)^2;
end
end

