function q_dot2=q_dot2_plan(prof,t)
% This function receives the method to calculate v and the time t 
% and calculat q_dot2 with the invers of the liniear jacobian
a = a_plan(prof,t);
q = q_plan(prof,t);
J = jacobian_mat(q);
J_l = J(1:3,1:3);
q_dot = q_dot_plan(prof,t);
% q_dot(1:2) = q_dot(1:2)*(pi/180);
J_dot = jacobian_mat_dot(q,q_dot);
J_dotl = J_dot(1:3,1:3);
q_dot2 = inv(J_l)*(a.'-J_dotl*q_dot(1:3).');
% q_dot2(1:2) = q_dot2(1:2)*(180/pi);
q_dot2 = q_dot2.';
end