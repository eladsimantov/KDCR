function J_dot = jacobian_mat_dot(q,q_dot)
% This function receives the positions and velocities of the joints,
% and returns the Jacobian derivative.
global H L 
l1 = 0;
l2 = 0;
theta1 = q(1);
theta2 = q(2);
d3 = q(3);
theta4 = 0;
theta5 = 0;
theta1_diff = q_dot(1);
theta2_diff = q_dot(2);
d3_diff = q_dot(3);
theta4_diff = 0;
theta5_diff = 0;
J_dot1 = [cos(theta1)*sin(theta2)*d3_diff-L*cos(theta1)*theta1_diff+...
    l1*cos(theta1)*cos(theta2)*theta2_diff+ l2*cos(theta1)*cos(theta2)*theta2_diff...
    - l1*sin(theta1)*sin(theta2)*theta1_diff - l2*sin(theta1)*sin(theta2)*theta1_diff...
    + cos(theta1)*cos(theta2)*d3*theta2_diff - sin(theta1)*sin(theta2)*d3*theta1_diff;...
    sin(theta1)*sin(theta2)*d3_diff - L*sin(theta1)*theta1_diff + ...
    l1*cos(theta1)*sin(theta2)*theta1_diff + l1*cos(theta2)*sin(theta1)*theta2_diff +...
    l2*cos(theta1)*sin(theta2)*theta1_diff + l2*cos(theta2)*sin(theta1)*theta2_diff...
    + cos(theta1)*sin(theta2)*d3*theta1_diff + cos(theta2)*sin(theta1)*d3*theta2_diff;...
    0;...
    0;...
    0;...
    0;];
J_dot2 = [cos(theta1)*cos(theta2)*theta1_diff*(l1 + l2 + d3) - ...
    sin(theta1)*(sin(theta2)*d3*theta2_diff - cos(theta2)*d3_diff + l1*sin(theta2)*theta2_diff + l2*sin(theta2)*theta2_diff);...
    cos(theta1)*(sin(theta2)*d3*theta2_diff - cos(theta2)*d3_diff +...
    l1*sin(theta2)*theta2_diff + l2*sin(theta2)*theta2_diff) + cos(theta2)*sin(theta1)*theta1_diff*(l1 + l2 + d3);...
    - sin(theta2)*d3_diff - cos(theta2)*d3*theta2_diff - l1*cos(theta2)*theta2_diff - l2*cos(theta2)*theta2_diff;...
    -sin(theta1)*theta1_diff;...
    cos(theta1)*theta1_diff;...
    0];
J_dot3 = [cos(theta1)*sin(theta2)*theta1_diff + cos(theta2)*sin(theta1)*theta2_diff;...
    sin(theta1)*sin(theta2)*theta1_diff - cos(theta1)*cos(theta2)*theta2_diff;...
    -sin(theta2)*theta2_diff;...
    0;...
    0;...
    0];
J_dot4 = [0;0;0;...
    cos(theta1)*sin(theta2)*theta1_diff + cos(theta2)*sin(theta1)*theta2_diff;...
    sin(theta1)*sin(theta2)*theta1_diff - cos(theta1)*cos(theta2)*theta2_diff;...
    -sin(theta2)*theta2_diff];
J_dot5 = [l2*sin(theta1)*theta1_diff;...
    -l2*cos(theta1)*theta1_diff;...
    0;...
    cos(theta1)*cos(theta2)*theta1_diff - sin(theta1)*sin(theta2)*theta2_diff;...
    cos(theta2)*sin(theta1)*theta1_diff + cos(theta1)*sin(theta2)*theta2_diff;...
    -cos(theta2)*theta2_diff];
% J_dot = [J_dot1 J_dot2 J_dot3 J_dot4 J_dot5];
J_dot = [J_dot1 J_dot2 J_dot3];
end