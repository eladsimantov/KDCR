function [] = Task_traj_plotter(profile_name, X, V, A, t)
%TRAJ_PLOTTER creates three figures of position velocity and acceleration
%for a given profile and its Trajectory data.
%   X - Waypoints with three columns and length t points.
%   V - Velocity. A - Acceleration.
% close all; 
figure('Name', 'Position Profile - '+ profile_name, 'NumberTitle', 'off');
subplot(3,1,1)
plot(t, X(:,1), LineStyle=":", LineWidth=1.5, Color="r")
grid on; title("Position over time - " + profile_name + " Velocity Profile"); 
ylabel('$x\ (m)$', 'Interpreter', 'latex');
subplot(3,1,2)
plot(t, X(:,2), LineStyle=":", LineWidth=1.5, Color="g")
grid on; ylabel('$y\ (m)$', 'Interpreter', 'latex');
subplot(3,1,3)
plot(t, X(:,3), LineStyle=":", LineWidth=1.5, Color="b")
grid on; ylabel('$z\ (m)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')
drawnow
pause(0.1)

figure('Name', 'Velocity Profile - '+ profile_name, 'NumberTitle', 'off');
subplot(3,1,1)
plot(t, V(:,1), LineStyle=":", LineWidth=1.5, Color="r")
grid on; title("Velocity over time - " + profile_name + " Velocity Profile"); 
ylabel('$\dot{x}\ (m/s)$', 'Interpreter', 'latex');
subplot(3,1,2)
plot(t, V(:,2), LineStyle=":", LineWidth=1.5, Color="g")
grid on; ylabel('$\dot{y}\ (m/s)$', 'Interpreter', 'latex');
subplot(3,1,3)
plot(t, V(:,3), LineStyle=":", LineWidth=1.5, Color="b")
grid on; ylabel('$\dot{z}\ (m/s)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')
drawnow
pause(0.1)

figure('Name', 'Acceleration Profile - '+ profile_name, 'NumberTitle', 'off');
subplot(3,1,1)
plot(t, A(:,1), LineStyle=":", LineWidth=1.5, Color="r")
grid on; title("Acceleration over time - " + profile_name + " Velocity Profile"); 
ylabel('$\ddot{x}\ (m/s^2)$', 'Interpreter', 'latex');
subplot(3,1,2)
plot(t, A(:,2), LineStyle=":", LineWidth=1.5, Color="g")
grid on; ylabel('$\ddot{y}\ (m/s^2)$', 'Interpreter', 'latex');
subplot(3,1,3)
plot(t, A(:,3), LineStyle=":", LineWidth=1.5, Color="b")
grid on; ylabel('$\ddot{z}\ (m/s^2)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')

end

