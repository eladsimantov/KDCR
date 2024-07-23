%% Test inv_kin and forward_kin
clear all; clc; close all;
global r L H
r=1; L=2; H=3;

%% Test FK and IK with specific task and find numerical errors
task = [2 1 deg2rad(45)]'; % [x; y; theta]
numerical_errors = test_fk_on_ik(task, 1e-12);

%% Draw all FK solutions for a given IK solution:
inv_kin_sols = inv_kin(task); % find the IK for the task variables
display(task)
display(inv_kin_sols)

% for all 8 possible IK sols, find the FK solutions to find the task
for j=1:8
    fprintf("Testing IK solution number: " + j)
    forward_kin_sol = forward_kin(inv_kin_sols(1:3,j));
    display(forward_kin_sol)
end


%% Draw all IK solutions for the given task:
for i=1:8
draw_parallel_robot(task, inv_kin_sols(1:3, i))
end

%% Singularity of IK

x_choice = 1;
sing_task_2 = [x_choice H-L-r/2*sqrt(3) 0]';
inv_kin_sols = inv_kin(sing_task_2); % 8 solutions, but really only 4 because two solution branches converged.
draw_parallel_robot(sing_task_2, inv_kin_sols(1:3,1))
draw_parallel_robot(sing_task_2, inv_kin_sols(1:3,3))
draw_parallel_robot(sing_task_2, inv_kin_sols(1:3,5))
draw_parallel_robot(sing_task_2, inv_kin_sols(1:3,7)) % skip every two because of converged solution branches.

% the other solutions types are not within the limits of 0<y<2..

%% Singularity of FK

% Solve numerically for x y theta with assumed ranges for y.
syms y1 d1 d2 d3 real
x1=1; theta=0;
assume(y1>0)
assume(y1<2)

% Eq 1,2,3 - take the three geometrical constraint equations.
% Eq 4 - take the determinant of Jx and equate to zero.
Eq1 = y1^2 + (d1 - x1)^2 - L^2 == 0;
Eq2 = (x1 - d2 + r*cos(theta))^2 + (y1 + r*sin(theta))^2 - L^2 == 0; 
Eq3 = L^2 == (y1 - H + r*sin(theta + pi/3))^2 + (x1 - d3 + r*cos(theta + pi/3))^2;
Eq4 = 4*d1*r*y1^2 + 4*d2*r*y1^2 - 8*d3*r*y1^2 + 4*H*r^2*y1 + 8*H*r*x1*y1 + 4*sqrt(3)*d1*r^2*y1 - 4*sqrt(3)*d3*r^2*y1 - 4*H*d1*r*y1 - 4*H*d2*r*y1 - 4*sqrt(3)*d1*d3*r*y1 + 4*sqrt(3)*d2*d3*r*y1 + 4*sqrt(3)*d1*r*x1*y1 - 4*sqrt(3)*d2*r*x1*y1==0;

digits(8)
sing_sol = vpasolve([Eq1, Eq2, Eq3, Eq4], [d1, d2, d3, y1]); % Solve 4 vars with 4 equations.
joints_sing_sols = [sing_sol.d1 sing_sol.d2 sing_sol.d3]; % grab joint values matrix.
task_sing_sols = [ones(6,1), sing_sol.y1, zeros(6,1)]; % grab task values matrix.

%% Tests of FK Singularities and Metrics.
% Test that these solutions are valid by nesting the IK on the task. 
for i=1:6
    task = task_sing_sols(i,1:3); % [x; y; theta] of one solution
    display(task)
    display(joints_sing_sols(i,1:3))
    inv_kin_sols = inv_kin(eval(task)); % find the IK for the task variables
    disp(double(inv_kin_sols))
end

% Test using Jx that these solutions bring det(Jx)=0.
display(joints_sing_sols)
display(task_sing_sols)
for i=1:6
    d_1 = joints_sing_sols(i,1);
    d_2 = joints_sing_sols(i,2);
    d_3 = joints_sing_sols(i,3);
    y_1 = task_sing_sols(i,2);
    det_Jx = @(d1, d2, d3, y1)(4*d1*r*y1^2 + 4*d2*r*y1^2 - 8*d3*r*y1^2 + 4*H*r^2*y1 + 8*H*r*x1*y1 + 4*sqrt(3)*d1*r^2*y1 - 4*sqrt(3)*d3*r^2*y1 - 4*H*d1*r*y1 - 4*H*d2*r*y1 - 4*sqrt(3)*d1*d3*r*y1 + 4*sqrt(3)*d2*d3*r*y1 + 4*sqrt(3)*d1*r*x1*y1 - 4*sqrt(3)*d2*r*x1*y1);
    Jx = @(d1,d2,d3,x1,y1)([2*d1 - 2*x1, -2*y1, 0; 2*d2 - 2*r - 2*x1, -2*y1, -2*r*y1; r - 2*d3 + 2*x1, 2*y1 - 2*H + sqrt(3)*r, r*y1 - H*r + sqrt(3)*d3*r - sqrt(3)*r*x1]);
    % det_singular = double(det_Jx(d_1, d_2, d_3, y_1))
    singular_direction = null(double(Jx(d_1, d_2, d_3, 1, y_1)));
    singular_direction = singular_direction ./ norm(singular_direction);
    % [vec,D] = eig(Jx(d_1, d_2, d_3, 1, y_1)) % another option.
    draw_parallel_robot_singularity_add_ins(double([1, y_1, 0]), double([d_1, d_2, d_3]), singular_direction)
    display(joints_sing_sols(i,1:3))
    display(singular_direction)
    display(double(Jx(d_1, d_2, d_3, 1, y_1))*singular_direction)
end


