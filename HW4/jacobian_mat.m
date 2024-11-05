function J = jacobian_mat(q)
%JACOBIAN_MAT Calculate the full Jacobian Matrix for given joint variable
%vector. 
%   q - joint position vector to compute full jacobian matrix.
global L
JL = [
    [q(3)*cos(q(1))*sin(q(2)) - L*sin(q(1)), q(3)*cos(q(2))*sin(q(1)), sin(q(1))*sin(q(2))];
    [L*cos(q(1)) + q(3)*sin(q(1))*sin(q(2)), -q(3)*cos(q(1))*cos(q(2)), -cos(q(1))*sin(q(2))];
    [0, -q(3)*sin(q(2)), cos(q(2))]];
JA = [
    [0, cos(q(1)), 0];
    [0, sin(q(1)), 0];
    [1, 0, 0]];
J = [JL; JA];
end


