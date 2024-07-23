function [X] = forward_kins(Q)
%FORWARD_KINS
%   Q – Matrix of row vectors of different joints parameters.
%   X – Matrix of position row vector of the tool.
[Nrows, Ncols] = size(Q);
X = zeros(Nrows, Ncols);
for row=1:Nrows
    q = Q(row, :);
    X(row, :) = forward_kin(q);
end
end

