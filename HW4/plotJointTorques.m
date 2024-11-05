function [] = plotJointTorques(plot_title, Tau, Tau_num, t)
figure('Name', plot_title, 'NumberTitle', 'off');
set(gcf, "Color", 'w'); 
subplot(3,1,1); hold on; box on;
plot(t, Tau(:,1), LineStyle="-", LineWidth=1.5)
plot(t, Tau_num(:,1), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; title(plot_title, "FontSize", 18); 
ylabel('${\tau}_1\ (Nm)$', 'Interpreter', 'latex'); hold off;
subplot(3,1,2); hold on; box on;
plot(t, Tau(:,2), LineStyle="-", LineWidth=1.5)
plot(t, Tau_num(:,2), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; ylabel('${\tau}_2\ (Nm)$', 'Interpreter', 'latex'); hold off;
subplot(3,1,3); hold on; box on;
plot(t, Tau(:,3), LineStyle="-", LineWidth=1.5)
plot(t, Tau_num(:,3), LineStyle=":", LineWidth=1.5)
legend(["Planned", "Dynamics"])
grid on; ylabel('$F_3\ (N)$', 'Interpreter', 'latex');
xlabel('$t\ (s)$', 'Interpreter', 'latex')

end

