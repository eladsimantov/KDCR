function q_dot=q_dot_plan(prof,t)
% This function receives a decision value for the velocity profile and the time t 
% and calculat q_dot with the invers of the liniear jacobian
v = v_plan(prof,t);
q = q_plan(prof,t);
J = jacobian_mat(q);
J_l = J(1:3,1:3);
q_dot = inv(J_l)*v.';
% q_dot(1:2) = q_dot(1:2)*(180/pi);
q_dot = q_dot.';
end