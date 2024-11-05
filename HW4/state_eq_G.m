function Xdot=state_eq_G(t,X)
% State simulation of the system
global Fe M_real
q1 = X(1);
q2 = X(3);
q3 = X(5);
q_dot_1 = X(2);
q_dot_2 = X(4);
q_dot_3 = X(6);
% Enforce joint 2 limits
if abs(q2) >= pi/2
    % Reach Hard Stop
   q2 = pi/2*sign(q2);
   % Check if velocity tries to increase q2 over the hard stop
   if sign(q2) == sign(q_dot_2)
       q_dot_2 = 0;
%         q_dot_2 = -q_dot_2; % Completely elastic collision
   end
end

% Calc the actual system dynamics (actual dynamics depend on real mass)
[H_mat,C_mat,G_mat] = dynamics_mat([q1,q2,q3],[q_dot_1,q_dot_2,q_dot_3], M_real);

tau = tau_cont([q1,q2,q3],[q_dot_1,q_dot_2,q_dot_3], t, "G");
J = jacobian_mat([q1 q2 q3]);
q_ddot = H_mat\(tau-C_mat*[q_dot_1;q_dot_2;q_dot_3]-G_mat+J'*Fe);
% q_ddot = H_mat\(-C_mat*[q_dot_1;q_dot_2;q_dot_3]-G_mat);
g1 = q_dot_1;
g2 = q_ddot(1);
g3 = q_dot_2;
g4 = q_ddot(2);
g5 = q_dot_3;
g6 = q_ddot(3);
Xdot = [g1;g2;g3;g4;g5;g6];
end