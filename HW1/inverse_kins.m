function [Q] = inverse_kins(X,elbows)
%INVERSE_KIN is a function to deal with a matrix of points and return a
%matrix of joint variables. It uses the simple inverse_kin function within.
%   X - Matrix of row position vectors of the tool.
%   elbows – Vector of decision values for the solutions in all times.
[Nrows, Ncols] = size(X);
Q = zeros(Nrows, Ncols);
for row=1:Nrows
    x = X(row, :);
    Q(row, :) = inverse_kin(x,elbows(row, :));
end
end

