function [Q_dot2] = q_dot2_plan(prof, t, elbows)
%Q_DOT2_PLAN Takes the profile and the time vector t and returns an
%acceleration matrix for all joints in all times t using the "Jacobian"
%method. It uses a_plan, jacobian_mat, jacobian_mat_dot functions within
%   t - time vector.
%   elbows - solution decision variables vector used in the IK solver.
%   q_dot2 – the second time derivative of the joints’ parameters vector.
%   Q_dot2 – the acceleration of joints' parameters Matrix of all times t.

Q_dot2 = zeros(length(t), 3); % init

Q = q_plan(prof, t, elbows); % grab the joint values matrix.
Q_dot = q_dot_plan(prof, t, elbows); % grab the joint vel matrix.
A = a_plan(prof, t); % grab the task accel matrix

for i=1:length(t)
    J = jacobian_mat(Q(i, :)); % grab jacobian at given joint position
    JL = J(1:3,1:3); % grab linear jacobian from full jacobian matrix
    J_dot = jacobian_mat_dot(Q(i, :), Q_dot(i, :)); % now get J_dot
    JL_dot = J_dot(1:3, 1:3); % grab only linear part.
    rhs = transpose(A(i,:)) - JL_dot*transpose(Q_dot(i,:)); % middleware
    Q_dot2(i, :) = transpose(JL\rhs);
end

end

