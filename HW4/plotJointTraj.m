function [] = plotJointTraj(plot_title, Q, Q_num, t)
%plotJointTraj Plot joint space position 
%   Q - Matrix of position of all three joints over time t (in RAD)
%   t - time vector
%   Q_num - data using numerical method 

figure('Name', plot_title, 'NumberTitle', 'off');
set(gcf, "Color", 'w'); 
subplot(3,1,1); hold on; box on;
plot(t, rad2deg(Q(:,1)), LineStyle="-", LineWidth=1.5)
plot(t, rad2deg(Q_num(:,1)), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; title(plot_title, "FontSize", 18); 
ylabel('${\theta}_1\ (deg)$', 'Interpreter', 'latex'); hold off;
subplot(3,1,2); hold on; box on;
plot(t, rad2deg(Q(:,2)), LineStyle="-", LineWidth=1.5)
plot(t, rad2deg(Q_num(:,2)), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; ylabel('${\theta}_2\ (deg)$', 'Interpreter', 'latex'); hold off;
subplot(3,1,3); hold on; box on;
plot(t, Q(:,3), LineStyle="-", LineWidth=1.5)
plot(t, Q_num(:,3), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; ylabel('$d_3\ (m)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')

end

