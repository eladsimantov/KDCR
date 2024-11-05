function J = jacobian_mat_simbolic(q)
% The jacobian_mat function calculates the full Jacobian matrix
% of the robot for specific joint values.
global H L l;
% J = zeros(6,3);
% Calulating J:
% The transformation matrixs:
theta1 = q(1,1);
theta2 = q(1,2);
d3 = q(1,3);
% theta4 = q(1,4);
% theta5 = q(1,5);
A01 = [...
    cos(theta1) -sin(theta1) 0 0;...
    sin(theta1) cos(theta1) 0 0;...
    0 0 1 H;...
    0 0 0 1];
A12 = [...
    1 0 0 L;...
    0 cos(theta2) -sin(theta2) 0;...
    0 sin(theta2) cos(theta2) 0;...
    0 0 0 1];
A23 = [...
    1 0 0 0;...
    0 1 0 0;...
    0 0 1 d3;...
    0 0 0 1];
% A34 = [...
%     cos(theta4) -sin(theta4) 0 0;...
%     sin(theta4) cos(theta4) 0 0;...
%     0 0 1 l1;...
%     0 0 0 1];
% A4t = [...
%     cos(theta5) 0 -sin(theta5) -sin(theta5)*l2;...
%     0 1 0 0;...
%     sin(theta5) 0 cos(theta5) cos(theta5)*l2;...
%     0 0 0 1];
% Calculating the transformation matrixs
A02 = A01*A12;
A0t = A02*A23;
% A04 = A03*A34;
% A0t = A04*A4t;
% calculating b and u:
p_t = A0t(1:3,4);
b0 = p_t;
b1 = p_t - A01(1:3,4);
% b2 = p_t - A02(1:3,4);
% b3 = p_t - A0t(1:3,4);
% b4 = p_t - A04(1:3,4);
u0 = [0;0;1];
u1 = A01(1:3,1:3)*[1;0;0];
u2 = A02(1:3,1:3)*[0;0;1];
% u3 = A03(1:3,1:3)*[0;0;1];
% u4 = A04(1:3,1:3)*[0;-1;0];
% J = [cross(u0,b0), cross(u1,b1), u2, cross(u3,b3), cross(u4,b4);...
%     u0, u1, [0;0;0], u3, u4];
J = [cross(u0,b0), cross(u1,b1), u2;...
    u0, u1, [0;0;0]];

end