function [H_mat,C_mat,G_mat]=dynamics_mat(q,q_dot)
global H L l
% calcuation mass
d = 0.015;
rou = 7800;
m1 = pi*(d/2)^2*H*rou;
m2 = pi*(d/2)^2*L*rou;
m3 = pi*(d/2)^2*l*rou;
M = 0.5;
g = 9.81;
theta_1 = q(1);
theta_2 = q(2);
d_3 = q(3);
q1 = q(1);
q2 = q(2);
q3 = q(3);
q_dot_1 = q_dot(1);
q_dot_2 = q_dot(2);
q_dot_3 = q_dot(3);
A01 = [...
    cos(theta_1) -sin(theta_1) 0 0;...
    sin(theta_1) cos(theta_1) 0 0;...
    0 0 1 H;...
    0 0 0 1];
A12 = [...
    1 0 0 L;...
    0 cos(theta_2) -sin(theta_2) 0;...
    0 sin(theta_2) cos(theta_2) 0;...
    0 0 0 1];
A23 = [...
    1 0 0 0;...
    0 1 0 0;...
    0 0 1 d_3;...
    0 0 0 1];


A02 = A01*A12;
A0t = A02*A23;

% Rotation matrixes
R01 = A01(1:3,1:3);
R02 = A02(1:3,1:3);
R0t = A0t(1:3,1:3);

% Origen location world frame
Ori_1 = A01(1:3,4);
Ori_2 = A02(1:3,4);
Ori_t = A0t(1:3,4);

% World frame vectors
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Center of mass in linck frame
r1 = H/2;
r2 = L/2;
r3 = l/2;
rc1 = r1*e3;
rc2 = r2*e1;
rc3 = (d_3-r3)*e3;

% Center of mass in world frame
Rc1 = rc1;
Rc2 = R01*rc2+Ori_1;
Rc3 = R02*rc3+Ori_2;
Rct = Ori_t;



% linear velocite
% Rc1_dot = jacobian(Rc1,q)*q_dot.';
% Rc2_dot =  jacobian(Rc2,q)*q_dot.';
% Rc3_dot =  jacobian(Rc3,q)*q_dot.';
% Rct_dot =  jacobian(Rct,q)*q_dot.';

% moment of inertia in link frame
I1 = (m1*H^2)/12;
I2 = (m2*L^2)/12;
I3 = (m3*l^2)/12;
% I_t = M*d_3^2;
I_tilde_1 = I1*[1,0,0;
    0,1,0;
    0,0,0];
I_tilde_2 = I2*[0,0,0;
    0,1,0;
    0,0,1];
I_tilde_3 = I3*[1,0,0;
    0,1,0;
    0,0,0];

% I_tilde_t = I_t*[1,0,0;
%     0,1,0;
%     0,0,0];
% moment of inertia in world frame
Inertia1 = I_tilde_1;
Inertia2 = R01*I_tilde_2*R01';
Inertia3 = R02*I_tilde_3*R02';
% Inertiat = R02*I_tilde_t*R02';
% % Angular velocity of the links
% omega1 = q_dot_1*e3;
% omega2 = omega1 + q_dot_2*(R01*e1);
% omega3 = omega2;

% Jacobian:
J = [q3*cos(q1)*sin(q2) - L*sin(q1), q3*cos(q2)*sin(q1), sin(q1)*sin(q2);...
    L*cos(q1) + q3*sin(q1)*sin(q2), -q3*cos(q1)*cos(q2), -cos(q1)*sin(q2);...
    0, cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2)), cos(q2);...
    0, cos(q1), 0;...
    0, sin(q1),0;...
    1, 0, 0];
J1 = [0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    1, 0, 0];
J2 = [-(L*sin(q1))/2,       0, 0;
    (L*cos(q1))/2,       0, 0;
    0,       0, 0;
    0, cos(q1), 0;
    0, sin(q1), 0;
    1,       0, 0];
J3 = [- L*sin(q1) - cos(q1)*sin(q2)*(l/2 - q3),-cos(q2)*sin(q1)*(l/2 - q3),  sin(q1)*sin(q2);
    L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3), cos(q1)*cos(q2)*(l/2 - q3), -cos(q1)*sin(q2);
    0, cos(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3)),          cos(q2);
    0,                                                                                             cos(q1),                0;
    0,                                                                                             sin(q1),                0;
    1,                                                                                                   0,                0];

% C_mat = [(M*q3^2*q_dot_2*sin(2*q2))/2 + (l^2*m3*q_dot_2*sin(2*q2))/6 + (m3*q3^2*q_dot_2*sin(2*q2))/2 + M*q3*q_dot_3*sin(q2)^2 - (l*m3*q_dot_3*sin(q2)^2)/2 + m3*q3*q_dot_3*sin(q2)^2 - (l*m3*q3*q_dot_2*sin(2*q2))/2, (M*q3^2*q_dot_1*sin(2*q2))/2 - L*M*q_dot_3*cos(q2) + (l^2*m3*q_dot_1*sin(2*q2))/6 + (m3*q3^2*q_dot_1*sin(2*q2))/2 - L*m3*q_dot_3*cos(q2) - (l*m3*q3*q_dot_1*sin(2*q2))/2 + L*M*q3*q_dot_2*sin(q2) - (L*l*m3*q_dot_2*sin(q2))/2 + L*m3*q3*q_dot_2*sin(q2), (l*m3*q_dot_1*(cos(q2)^2 - 1))/2 - M*q3*q_dot_1*(cos(q2)^2 - 1) - L*M*q_dot_2*cos(q2) - m3*q3*q_dot_1*(cos(q2)^2 - 1) - L*m3*q_dot_2*cos(q2);
%  -(q_dot_1*sin(2*q2)*(3*M*q3^2 + l^2*m3 + 3*m3*q3^2 - 3*l*m3*q3))/6,(q_dot_3*(2*M*q3 - l*m3 + 2*m3*q3))/2,(q_dot_2*(2*M*q3 - l*m3 + 2*m3*q3))/2;
% -(q_dot_1*sin(q2)^2*(2*M*q3 - l*m3 + 2*m3*q3))/2,-(q_dot_2*(2*M*q3 - l*m3 + 2*m3*q3))/2,0];
C_mat = [- q_dot_3*(M*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + m3*cos(q1)*sin(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - m3*sin(q1)*sin(q2)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) - q_dot_2*(M*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) - m3*cos(q1)*cos(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3) + m3*cos(q2)*sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3)),- q_dot_2*(M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q1)^2*cos(q2)*sin(q2) + q3*cos(q2)*sin(q1)^2*sin(q2)) - M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + cos(q2)*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2) - q3*cos(q1)^2*cos(q2)^2 + q3*cos(q1)^2*sin(q2)^2 - q3*cos(q2)^2*sin(q1)^2 + q3*sin(q1)^2*sin(q2)^2) - M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(q3*cos(q1)^2*cos(q2)^2 - sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q2)^2*sin(q1)^2) - M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) - sin(q2)*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2) + 2*q3*cos(q1)^2*cos(q2)*sin(q2) + 2*q3*cos(q2)*sin(q1)^2*sin(q2)) + (L^2*m2*cos(q1)*sin(q1))/12 + m3*cos(q1)*sin(q2)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3) + m3*sin(q1)*sin(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3)) - q_dot_3*((M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2)))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2))/2 - (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2))/2 + (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)))/2 + m3*cos(q1)*cos(q2)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3)) + m3*cos(q2)*sin(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3))) - q_dot_1*(M*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) - m3*cos(q1)*cos(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3) + m3*cos(q2)*sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3)),- q_dot_1*(M*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + m3*cos(q1)*sin(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - m3*sin(q1)*sin(q2)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) - q_dot_2*((M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2)))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2))/2 - (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2))/2 + (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)))/2 + m3*cos(q1)*cos(q2)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3)) + m3*cos(q2)*sin(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)));...
    q_dot_1*(M*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) - m3*cos(q1)*cos(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3) + m3*cos(q2)*sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))*(l/2 - q3)) - q_dot_3*((M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2))/2 + (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)))/2) + (L^2*m2*q_dot_2*cos(q1)*sin(q1))/12,(L^2*m2*q_dot_1*cos(q1)*sin(q1))/12 - q_dot_2*(M*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q1)^2*cos(q2)*sin(q2) + q3*cos(q2)*sin(q1)^2*sin(q2))*(sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + cos(q2)*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2) - q3*cos(q1)^2*cos(q2)^2 + q3*cos(q1)^2*sin(q2)^2 - q3*cos(q2)^2*sin(q1)^2 + q3*sin(q1)^2*sin(q2)^2) - m3*(cos(q1)^2*cos(q2)*(l/2 - q3) + cos(q2)*sin(q1)^2*(l/2 - q3))*(cos(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) + M*(q3*cos(q1)^2*cos(q2)^2 - sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q2)^2*sin(q1)^2)*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) - sin(q2)*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2) + 2*q3*cos(q1)^2*cos(q2)*sin(q2) + 2*q3*cos(q2)*sin(q1)^2*sin(q2)) + m3*cos(q1)^2*cos(q2)*sin(q2)*(l/2 - q3)^2 + m3*cos(q2)*sin(q1)^2*sin(q2)*(l/2 - q3)^2) - q_dot_3*(m3*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) - M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q1)^2*cos(q2)*sin(q2) + q3*cos(q2)*sin(q1)^2*sin(q2)) - M*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2)*(q3*cos(q1)^2*cos(q2)^2 - sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q2)^2*sin(q1)^2) + (m3*cos(q1)^2*cos(q2)^2*(l - 2*q3))/2 + (m3*cos(q2)^2*sin(q1)^2*(l - 2*q3))/2),q_dot_3*(M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2) - M*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2) + M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2) + m3*cos(q2)*sin(q2) + M*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2)*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)) - m3*cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)) - q_dot_2*(m3*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) - M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q1)^2*cos(q2)*sin(q2) + q3*cos(q2)*sin(q1)^2*sin(q2)) - M*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2)*(q3*cos(q1)^2*cos(q2)^2 - sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q2)^2*sin(q1)^2) + (m3*cos(q1)^2*cos(q2)^2*(l - 2*q3))/2 + (m3*cos(q2)^2*sin(q1)^2*(l - 2*q3))/2) - q_dot_1*((M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2))/2 + (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)))/2);...
    q_dot_1*(M*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + m3*cos(q1)*sin(q2)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - m3*sin(q1)*sin(q2)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) + q_dot_2*((M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2))/2 + (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)))/2),q_dot_1*((M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2)))/2 - (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2))/2 + (M*(cos(q1)*cos(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + cos(q2)*sin(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2))/2 + (M*(cos(q1)*sin(q2)*(L*cos(q1) + q3*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(L*sin(q1) - q3*cos(q1)*sin(q2)))*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2)))/2) - q_dot_2*(M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) - sin(q2)*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2) + 2*q3*cos(q1)^2*cos(q2)*sin(q2) + 2*q3*cos(q2)*sin(q1)^2*sin(q2)) - m3*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2)*(cos(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) + M*(q3*cos(q1)^2*cos(q2)^2 - sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q2)^2*sin(q1)^2)*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2) + M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q2)*sin(q1)^2*sin(q2))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q1)^2*cos(q2)*sin(q2) + q3*cos(q2)*sin(q1)^2*sin(q2)) + M*(sin(q2)*(sin(q2)*cos(q1)^2 + sin(q2)*sin(q1)^2) + cos(q1)^2*cos(q2)^2 + cos(q2)^2*sin(q1)^2)*(q3*cos(q1)^2*cos(q2)^2 - sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q2)^2*sin(q1)^2) + m3*sin(q2)*(cos(q1)*(L*sin(q1) + cos(q1)*sin(q2)*(l/2 - q3)) - sin(q1)*(L*cos(q1) - sin(q1)*sin(q2)*(l/2 - q3))) - M*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2))*(cos(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + q3*cos(q1)^2*cos(q2)*sin(q2) + q3*cos(q2)*sin(q1)^2*sin(q2)) - m3*cos(q2)*(cos(q1)^2*cos(q2)*(l/2 - q3) + cos(q2)*sin(q1)^2*(l/2 - q3)) + M*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2)*(sin(q2)*(cos(q1)*(L*sin(q1) - q3*cos(q1)*sin(q2)) - sin(q1)*(L*cos(q1) + q3*sin(q1)*sin(q2))) + cos(q2)*(q3*cos(q2)*cos(q1)^2 + q3*cos(q2)*sin(q1)^2) - q3*cos(q1)^2*cos(q2)^2 + q3*cos(q1)^2*sin(q2)^2 - q3*cos(q2)^2*sin(q1)^2 + q3*sin(q1)^2*sin(q2)^2) + m3*cos(q1)^2*cos(q2)^2*(l/2 - q3) - m3*cos(q1)^2*sin(q2)^2*(l/2 - q3) + m3*cos(q2)^2*sin(q1)^2*(l/2 - q3) - m3*sin(q1)^2*sin(q2)^2*(l/2 - q3) - (m3*cos(q1)^2*cos(q2)^2*(l - 2*q3))/2 - (m3*cos(q2)^2*sin(q1)^2*(l - 2*q3))/2) + q_dot_3*(M*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2) - M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2) - m3*cos(q2)*sin(q2) + m3*cos(q1)^2*cos(q2)*sin(q2) + m3*cos(q2)*sin(q1)^2*sin(q2)),q_dot_2*(M*(2*cos(q1)^2*cos(q2)*sin(q2) - 2*cos(q2)*sin(q2) + 2*cos(q2)*sin(q1)^2*sin(q2))*(cos(q1)^2*sin(q2)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2) - M*(cos(q1)^2*cos(q2)*sin(q2) - cos(q2)*sin(q2) + cos(q2)*sin(q1)^2*sin(q2))*(- cos(q1)^2*cos(q2)^2 + cos(q1)^2*sin(q2)^2 - cos(q2)^2*sin(q1)^2 + cos(q2)^2 + sin(q1)^2*sin(q2)^2 - sin(q2)^2) - m3*cos(q2)*sin(q2) + m3*cos(q1)^2*cos(q2)*sin(q2) + m3*cos(q2)*sin(q1)^2*sin(q2))];

H_mat = m1*J1(1:3,:)'*J1(1:3,:)+m2*J2(1:3,:)'*J2(1:3,:)+(m3)*J3(1:3,:)'*J3(1:3,:) + M*J(1:3,:)'*J(1:3,:)+ J1(4:6,:)'*Inertia1*J1(4:6,:) + J2(4:6,:)'*Inertia2*J2(4:6,:)+J3(4:6,:)'*Inertia3*J3(4:6,:);
G_mat = g*(m1*J1(1:3,:)'+m2*J2(1:3,:)'+(m3)*J3(1:3,:)'+M*J(1:3,:)')*e3;
  
end