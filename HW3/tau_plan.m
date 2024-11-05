function tau=tau_plan(prof,t)
global H L l
q=q_plan(prof,t);
q_dot=q_dot_plan(prof,t);
q_dot2=q_dot2_plan(prof,t);
[H_mat,C_mat,G_mat] = dynamics_mat(q,q_dot);
tau = H_mat*q_dot2'+C_mat*q_dot'+G_mat;
end