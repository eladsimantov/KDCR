%% ----------------- Main -------------------- %
clc; close all; clearvars;
% Set global variables
global H L l elbow1 elbow2 elbow3 Fe M_real M_control
M_control = 0.5; 
M_real = 0.5;
H           =       0.20;
L           =       0.10;
l           =       0.60;
elbow1      =       1;
elbow2      =       1;
elbow3      =       1;
Fe          =       zeros(6,1);

% Default plot options
set(0, 'DefaultAxesFontSize', 15);  
set(0, 'DefaultLegendFontSize', 14); 
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0, 'DefaultAxesFontName', 'LaTeX');       

% Prep simulation
dt = 0.001;
tspan = 0:dt:3;
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-14, 'MaxStep',dt); % found to yield best results, to avoid oscillations due to step size and interpolations in ode solver.

% Get the planned x, q and tau
q_planned = reshape(cell2mat(arrayfun(@(t)(q_plan(3,t)), tspan,"UniformOutput",false)), 3, []).';
q_dot_planned = reshape(cell2mat(arrayfun(@(t)(q_dot_plan(3,t)), tspan,"UniformOutput",false)), 3, []).';
x_planned = reshape(cell2mat(arrayfun(@(t)(x_plan(3,t)), tspan,"UniformOutput",false)), 3, []).';
tau_planned = reshape(cell2mat(arrayfun(@(t)(tau_plan(3,t)), tspan,"UniformOutput",false)), 3, []).';

% Add some error to the IC to see the error response
x_0_err = [0, 0, 0.01]; % error in initial position of 1cm upward +z_0
x_0 = x_plan(3,0); 
q_0 = inverse_kin(x_0 + x_0_err,[elbow1 elbow2 elbow3],1); % includes the error
X0 = [q_0(1);0;q_0(2);0;q_0(3);0]; % State vector ICs (with the error in z axis)

%% ------------------------ Inverse Dynamics ---------------------------- %
clc; close all;

% run simulation
[t,X_num] = ode45(@(t_sim, X_sim) state_eq_ID(t_sim, X_sim),tspan,X0,options); % Inverse Dynamics

% Get the numerical x, q, qdot and tau from simulation
q_num = [X_num(:,1) X_num(:,3) X_num(:,5)];
qdot_num = [X_num(:,2) X_num(:,4) X_num(:,6)];
tau_num = zeros(size(q_num));
for i=1:length(t)
    tau_num(i, :) = tau_cont(q_num(i, :), qdot_num(i,:), t(i), "ID");
end
x_num = forward_kin([X_num(:,1) X_num(:,3) X_num(:,5)]);

% Calc the errors in x, q and tau
eq = q_num - q_planned;             % numerical q minus desired q
ex = x_num - x_planned;             % numerical x minus desired x
etau = tau_num - tau_planned;       % numerical tau minus desired tau

% Plots
plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off'); subplot(3,1,1); ylim([-1e-1, 1e-1]);
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)
plotTaskTraj("Task Values", x_planned, x_num, t)
plotTaskTraj("Task Values Error", zeros(size(ex)), ex, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');

%% ID - Dynamic Simulation (no load mass) 
clc; close all;
M_real = 0;
% Run simulation for ID again with M=0 for Q3
[t,X_num] = ode45(@(t_sim, X_sim) state_eq_ID(t_sim, X_sim),tspan,X0,options); % Inverse Dynamics

% Get the numerical x, q, qdot and tau from simulation
q_num = [X_num(:,1) X_num(:,3) X_num(:,5)];
qdot_num = [X_num(:,2) X_num(:,4) X_num(:,6)];
tau_num = zeros(size(q_num));
for i=1:length(t)
    tau_num(i, :) = tau_cont(q_num(i, :), qdot_num(i,:), t(i), "ID");
end
x_num = forward_kin([X_num(:,1) X_num(:,3) X_num(:,5)]);

% Calc the errors in x, q and tau
eq = q_num - q_planned;             % numerical q minus desired q
ex = x_num - x_planned;             % numerical x minus desired x
etau = tau_num - tau_planned;       % numerical tau minus desired tau

% Plots
plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)
plotTaskTraj("Task Values", x_planned, x_num, t)
plotTaskTraj("Task Values Error", zeros(size(ex)), ex, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');

%% ------------------------ PD + G ---------------------------- %
clc; close all;

% run simulation
M_real = 0.5;
[t,X_num] = ode45(@(t_sim, X_sim) state_eq_G(t_sim, X_sim),tspan,X0,options); 

% Get the numerical x, q, qdot and tau from simulation
q_num = [X_num(:,1) X_num(:,3) X_num(:,5)];
qdot_num = [X_num(:,2) X_num(:,4) X_num(:,6)];
tau_num = zeros(size(q_num));
for i=1:length(t)
    tau_num(i, :) = tau_cont(q_num(i, :), qdot_num(i,:), t(i), "G");
end
x_num = forward_kin([X_num(:,1) X_num(:,3) X_num(:,5)]);

% Calc the errors in x, q and tau
eq = q_num - q_planned;             % numerical q minus desired q
ex = x_num - x_planned;             % numerical x minus desired x
etau = tau_num - tau_planned;       % numerical tau minus desired tau

% Plots
plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off'); 
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)
plotTaskTraj("Task Values", x_planned, x_num, t)
plotTaskTraj("Task Values Error", zeros(size(ex)), ex, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');

%% PD+G - Dynamic simulation (no load mass)
clc; close all;

% Run simulation for ID again with M=0 for Q3
M_real = 0;
[t,X_num] = ode45(@(t_sim, X_sim) state_eq_G(t_sim, X_sim),tspan,X0,options); % Inverse Dynamics

% Get the numerical x, q, qdot and tau from simulation
q_num = [X_num(:,1) X_num(:,3) X_num(:,5)];
qdot_num = [X_num(:,2) X_num(:,4) X_num(:,6)];
tau_num = zeros(size(q_num));
for i=1:length(t)
    tau_num(i, :) = tau_cont(q_num(i, :), qdot_num(i,:), t(i), "ID");
end
x_num = forward_kin([X_num(:,1) X_num(:,3) X_num(:,5)]);

% Calc the errors in x, q and tau
eq = q_num - q_planned;             % numerical q minus desired q
ex = x_num - x_planned;             % numerical x minus desired x
etau = tau_num - tau_planned;       % numerical tau minus desired tau

% Plots
plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)


%% ------------------------ PID ---------------------------- %
clc; close all;
% Get a real mass between 0 and 1.
rng(1);
M_real = rand;
global KpPID KiPID KdPID
KpPID = [1300 0 0; 0 1700 0; 0 0 1100];
KdPID = [50 0 0; 0 60 0; 0 0 50];
KiPID = [100 0 0; 0 4000 0; 0 0 2000];

% Create new augmented state vector
% X_PID = [∫e1; ∫e2; ∫e3; e1; e2; e3; edot1; edot2; edot3]
q0_plan = q_plan(3, 0);
q0dot_plan = q_dot_plan(3, 0);
X0_PID = [zeros(3,1); X0(1)-q0_plan(1); X0(3)-q0_plan(2); X0(5)-q0_plan(3); X0(2)-q0dot_plan(1); X0(4)-q0dot_plan(2); X0(6)-q0dot_plan(3)];

% run simulation (with augmented state vector !!!)
global X_numPID % NO CHOICE BUT TO SAVE THE VALUES FOR TAU_CONT TO USE LATER...
% YES, YES, Think about it.. If you use a state equation in ode45 and you
% wish to calculate tau in each step post/irrespectively of simulation, so there is no choice but to save the simulated results.
% one does not simply calculate tau using t because ode45 will do it according to the RK numerical method. 
% The simulation will not use your "t" timesteps. "persistent" won't do the
% job correctly, it will only let you approximate tau, but not calculate
% the exact values used in the state equation.
[t,X_numPID] = ode45(@(t_sim, X_sim) state_eq_PID(t_sim, X_sim),tspan,X0_PID,options); 

% Grab data
eq = X_numPID(:, 4:6);
eqdot = X_numPID(:, 7:9);
q_num = eq + q_planned;
qdot_num = eqdot + q_dot_planned;
x_num = forward_kin(q_num);
ex = x_num - x_planned;
tau_num = zeros(size(q_num));
for i=1:length(t)
    tau_num(i,:) = -KpPID*X_numPID(i, 4:6).' - KdPID*X_numPID(i, 7:9).' - KiPID*X_numPID(i, 1:3).';
end
etau = tau_num - tau_planned;

plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off'); 
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)
plotTaskTraj("Task Values", x_planned, x_num, t)
plotTaskTraj("Task Values Error", zeros(size(ex)), ex, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');


%% PID dynamic simulation (no load mass)
clc; close all;
M_real = 0;
[t,X_num] = ode45(@(t_sim, X_sim) state_eq_PID(t_sim, X_sim),tspan,X0_PID,options); 
eq = X_num(:, 4:6);
eqdot = X_num(:, 7:9);
q_num = eq + q_planned;
qdot_num = eqdot + q_dot_planned;
x_num = forward_kin(q_num);
ex = x_num - x_planned;
tau_num = zeros(size(q_num));
for i=1:length(t)
    tau_num(i,:) = -KpPID*X_num(i, 4:6).' - KdPID*X_num(i, 7:9).' - KiPID*X_num(i, 1:3).';
end
etau = tau_num - tau_planned;

plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off'); 
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)



%% ------------------------ MINMAX ---------------------------- %
clc; close all;
rng(2)
M_real = rand;
[t,X_num] = ode45(@(t_sim, X_sim) state_eq_MINMAX(t_sim, X_sim),tspan,X0,options); 

% Get the numerical x, q, qdot and tau from simulation
q_num = [X_num(:,1) X_num(:,3) X_num(:,5)];
qdot_num = [X_num(:,2) X_num(:,4) X_num(:,6)];
tau_num = zeros(size(q_num));

for i=1:length(t)
    tau_num(i, :) = tau_cont(q_num(i, :), qdot_num(i,:), t(i), "MINMAX");
end
x_num = forward_kin([X_num(:,1) X_num(:,3) X_num(:,5)]);

% Calc the errors in x, q and tau
eq = q_num - q_planned;             % numerical q minus desired q
ex = x_num - x_planned;             % numerical x minus desired x
etau = tau_num - tau_planned;       % numerical tau minus desired tau

% Plots
plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off'); 
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)
plotTaskTraj("Task Values", x_planned, x_num, t)
plotTaskTraj("Task Values Error", zeros(size(ex)), ex, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');


%% MINMAX dynamic simulation (no load mass)
clc; close all;
M_real = 0;

[t,X_num] = ode45(@(t_sim, X_sim) state_eq_MINMAX(t_sim, X_sim),tspan,X0,options); 

% Get the numerical x, q, qdot and tau from simulation
q_num = [X_num(:,1) X_num(:,3) X_num(:,5)];
qdot_num = [X_num(:,2) X_num(:,4) X_num(:,6)];
tau_num = zeros(size(q_num));

for i=1:length(t)
    tau_num(i, :) = tau_cont(q_num(i, :), qdot_num(i,:), t(i), "MINMAX");
end
x_num = forward_kin([X_num(:,1) X_num(:,3) X_num(:,5)]);

% Calc the errors in x, q and tau
eq = q_num - q_planned;             % numerical q minus desired q
ex = x_num - x_planned;             % numerical x minus desired x
etau = tau_num - tau_planned;       % numerical tau minus desired tau

% Plots
plotJointTraj("Joint Values", q_planned, q_num, t)
plotJointTraj("Joint Values Error", zeros(size(eq)), eq, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off'); 
plotJointTorques("Joint Torques and Forces", tau_planned, tau_num, t)
plotJointTorques("Joint Torques and Forces Error", zeros(size(etau)), etau, t); set(findall(gcf, 'Type', 'legend'), 'Visible', 'off');
plotErrorNorm(x_num, x_planned, t)

