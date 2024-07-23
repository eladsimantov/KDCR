function J_dot = jacobian_mat_dot(q, q_dot)
%JACOBIAN_MAT Recieves joint position and velocity vector and returns the
%Jacobian time derivative.
%   J_dot – the jacobian time derivative. 
%   q - the joint current positions.
%   q_dot – the joint current velocities.
global L
J_dot = [
    [cos(q(1))*sin(q(2))*q_dot(3) - L*cos(q(1))*q_dot(1) + cos(q(1))*cos(q(2))*q(3)*q_dot(2) - sin(q(1))*sin(q(2))*q(3)*q_dot(1), cos(q(2))*sin(q(1))*q_dot(3) + cos(q(1))*cos(q(2))*q(3)*q_dot(1) - sin(q(1))*sin(q(2))*q(3)*q_dot(2), cos(q(1))*sin(q(2))*q_dot(1) + cos(q(2))*sin(q(1))*q_dot(2)];
    [sin(q(1))*sin(q(2))*q_dot(3) - L*sin(q(1))*q_dot(1) + cos(q(1))*sin(q(2))*q(3)*q_dot(1) + cos(q(2))*sin(q(1))*q(3)*q_dot(2), cos(q(2))*sin(q(1))*q(3)*q_dot(1) - cos(q(1))*cos(q(2))*q_dot(3) + cos(q(1))*sin(q(2))*q(3)*q_dot(2), sin(q(1))*sin(q(2))*q_dot(1) - cos(q(1))*cos(q(2))*q_dot(2)];
    [0, -sin(q(2))*q_dot(3)-cos(q(2))*q(3)*q_dot(2), -sin(q(2))*q_dot(2)]
    [0, -sin(q(1))*q_dot(1), 0];
    [0, cos(q(1))*q_dot(1), 0];
    [0, 0, 0]]; 
end

