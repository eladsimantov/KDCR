function [] = plotTaskTraj(plot_name, Xdes, X, t)
%TASK_TRAJ_PLOTTER Plot joint space position 
%   t - time vector

figure('Name', plot_name, 'NumberTitle', 'off');
set(gcf, "Color", 'w'); 

ax1 = subplot(3,1,1); hold on; box on;
plot(t, rad2deg(Xdes(:,1)), LineStyle="-", LineWidth=1.5)
plot(t, rad2deg(X(:,1)), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; title(plot_name, "FontSize", 18); 
ylabel('$x\ (m)$', 'Interpreter', 'latex'); hold off;

ax2 = subplot(3,1,2); hold on; box on;
plot(t, rad2deg(Xdes(:,2)), LineStyle="-", LineWidth=1.5)
plot(t, rad2deg(X(:,2)), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; ylabel('$y\ (m)$', 'Interpreter', 'latex'); hold off;

ax3 = subplot(3,1,3); hold on; box on;
plot(t, Xdes(:,3), LineStyle="-", LineWidth=1.5)
plot(t, X(:,3), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; ylabel('$z\ (m)$', 'Interpreter', 'latex');

linkaxes([ax1 ax2 ax3],'x')
xlabel('$t\ (s)$', 'Interpreter', 'latex'); 
end

