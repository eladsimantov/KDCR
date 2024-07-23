function [Q_dot] = q_dot_plan(prof, t, elbows)
%Q_DOT_PLAN takes a profile and a time vector and returns the joint 
%velocity using the "Jacobian" method. It uses the v_plan function along
%with the Jacobian_mat function.
%   t - the time vector, where i=length(t).
%   elbows - solution decision variables vector used in the IK solver.
%   q_dot – the joints velocity at time t(i) 
%   Q_dot - returned value - the joint velocity Matrix for all times t.

Q_dot = zeros(length(t), 3); % init

Q = q_plan(prof, t, elbows); % grab the joint values matrix
V = v_plan(prof, t); % grab the task velocity matrix

for i=1:length(t)
    J = jacobian_mat(Q(i, :)); % grab jacobian at given joint position
    JL = J(1:3,1:3); % grab linear jacobian from full jacobian matrix
    Q_dot(i, :) = transpose(JL\transpose(V(i, :))); % insert the joint velocity from equation into Qdot matrix at i-th row

end

end

    
    