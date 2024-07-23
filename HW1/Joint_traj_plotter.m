function [] = Joint_traj_plotter(profile_name, Q, Qdot_num, Qdot_jac, Qddot_num, Qddot_jac, t)
%JOINT_TRAJ_PLOTTER Plot joint space position velocity and acceleration in
%three figures. velocity and acceleration with two methods - one in 
%numerical method and one by using jacobian relations and the chain rule.
%   Q - Matrix of position of all three joints over time t (in RAD)
%   Qdot - Matrix of joint velocity (in RAD/s)
%   Qddot - Matrix of joint accel (in RAD/s2)
%   t - time vector
%   num - data using numerical method 
%   jac - data using the jacobian relation and chain rule


% close all; 
figure('Name', 'Joint Positions - '+ profile_name, 'NumberTitle', 'off');
subplot(3,1,1)
plot(t, rad2deg(Q(:,1)), LineStyle=":", LineWidth=1.5)
grid on; title("Joint Positions - " + profile_name + " Velocity Profile"); 
ylabel('$\theta_1\ (deg)$', 'Interpreter', 'latex');
subplot(3,1,2)
plot(t, rad2deg(Q(:,2)), LineStyle=":", LineWidth=1.5)
grid on; ylabel('$\theta_2\ (deg)$', 'Interpreter', 'latex');
subplot(3,1,3)
plot(t, Q(:,3), LineStyle=":", LineWidth=1.5)
grid on; ylabel('$d_3\ (m)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')
drawnow
pause(0.1)

figure('Name', 'Joint Velocities - '+profile_name, 'NumberTitle', 'off');
subplot(3,1,1)
plot(t, [rad2deg(Qdot_jac(:,1)), rad2deg(Qdot_num(:,1))], LineStyle=":", LineWidth=1.5)
legend(["Jacobians", "Numeric"])
grid on; title("Joint Velocities - " + profile_name + " Velocity Profile"); 
ylabel('$\dot{\theta}_1\ (deg/s)$', 'Interpreter', 'latex');
subplot(3,1,2)
plot(t, [rad2deg(Qdot_jac(:,2)), rad2deg(Qdot_num(:,2))], LineStyle=":", LineWidth=1.5)
legend(["Jacobians", "Numeric"])
grid on; ylabel('$\dot{\theta}_2\ (deg/s)$', 'Interpreter', 'latex');
subplot(3,1,3)
plot(t, [Qdot_jac(:,3), Qdot_num(:,3)], LineStyle=":", LineWidth=1.5)
legend(["Jacobians", "Numeric"])
grid on; ylabel('$\dot{d}_3\ (m/s)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')
drawnow
pause(0.1)


figure('Name', 'Joint Acceleration - '+profile_name, 'NumberTitle', 'off');
subplot(3,1,1)
plot(t, [rad2deg(Qddot_jac(:,1)), rad2deg(Qddot_num(:,1))], LineStyle=":", LineWidth=1.5)
legend(["Jacobians", "Numeric"])
grid on; title("Joint Acceleration - " + profile_name + " Velocity Profile"); 
ylabel('$\ddot{\theta}_1\ (deg/s^2)$', 'Interpreter', 'latex');
subplot(3,1,2)
plot(t, [rad2deg(Qddot_jac(:,2)), rad2deg(Qddot_num(:,2))], LineStyle=":", LineWidth=1.5)
legend(["Jacobians", "Numeric"])
grid on; ylabel('$\ddot{\theta}_2\ (deg/s^2)$', 'Interpreter', 'latex');
subplot(3,1,3)
plot(t, [Qddot_jac(:,3), Qddot_num(:,3)], LineStyle=":", LineWidth=1.5)
legend(["Jacobians", "Numeric"])
grid on; ylabel('$\ddot{d}_3\ (m/s^2)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')





end

