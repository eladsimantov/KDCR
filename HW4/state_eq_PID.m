function Xdot=state_eq_PID(t,X)
% State simulation of the system
global Fe M_real KpPID KiPID KdPID
% global eI
% if isempty(eI) || t==0
%     eI=[0;0;0];
% end

eI = X(1:3);
e = X(4:6);
edot = X(7:9);
dt=0.001;

q_des = q_plan(3, t);
qdot_des = q_dot_plan(3, t);
qdot2_des = q_dot2_plan(3,t);

q = e' + q_des;
qdot = edot' + qdot_des;
q1 = q(1);
q2 = q(2);
q3 = q(3);
q_dot_1 = qdot(1);
q_dot_2 = qdot(2);
q_dot_3 = qdot(3);

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

% Calc the system dynamics after checking q2 limits
[H_mat,C_mat,G_mat] = dynamics_mat([q1,q2,q3],[q_dot_1,q_dot_2,q_dot_3], M_real);
tau =   -KpPID*e - KdPID*edot - KiPID*eI;
% tau = tau_cont([q1,q2,q3],[q_dot_1,q_dot_2,q_dot_3], t, "PID");
% disp(tau - tau_c)

J = jacobian_mat([q1 q2 q3]);
q_ddot = H_mat\(tau-C_mat*[q_dot_1;q_dot_2;q_dot_3]-G_mat+J'*Fe);
eddot = q_ddot - qdot2_des';
Xdot = [e; edot; eddot];
end