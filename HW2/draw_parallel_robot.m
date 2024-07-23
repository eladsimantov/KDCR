function [] = draw_parallel_robot(x,q)
%DRAW_PARALLEL_ROBOT accepts both task and joint variables, plots the robot
%and saves the figure

global r L H

% Define Points
d1 = q(1); d2 = q(2); d3 = q(3);
x1 = x(1); y1 = x(2); theta = x(3);
x2 = x1 + r*cos(theta);
y2 = y1 + r*sin(theta);
x3 = x1 + r*cos(theta+pi/3);
y3 = y1 + r*sin(theta+pi/3);

figure()
hold on

% --------------------- %
% Draw Prizmatic links, some overlap so make them transparent & rectangles
if d1 >= 0
    rectangle('Position', [0, -0.04, d1, 0.08], 'FaceColor', [1, 0, 0, 0.1]);
else 
    rectangle('Position', [d1, -0.04, -d1, 0.08], 'FaceColor', [1, 0, 0, 0.1]);
end
rectangle('Position', [0, -0.02, d2, 0.04], 'FaceColor', [0, 0, 1, 0.1]);
if d3 >= 0 
    rectangle('Position', [0, H-0.02, d3, 0.04], 'FaceColor', [0, 1, 0, 0.5]);
else
    rectangle('Position', [d3, H-0.02, -d3, 0.04], 'FaceColor', [0, 1, 0, 0.5]);
end


% Draw Prizmatic joints as rectangles + circles - to show sliding mechanism 
width = 0.35; height = 0.25;
radius = height/5;
rectangle('Position', [d1-width/2, -height/2, width, height], 'EdgeColor', 'k', 'LineWidth', 1, FaceColor=[1,1,1,1]); 
rectangle('Position', [d2-width/2, -height/2, width, height], 'EdgeColor', 'k', 'LineWidth', 1, FaceColor=[1,1,1,1]); 
rectangle('Position', [d3-width/2, H-height/2, width, height], 'EdgeColor', 'k', 'LineWidth', 1, FaceColor=[1,1,1,1]);
rectangle('Position', [d1-radius, -radius, 2*radius, 2*radius], 'Curvature', [1, 1], EdgeColor='k',LineWidth=1, FaceColor=[0,0,0,1]);
rectangle('Position', [d2-radius, -radius, 2*radius, 2*radius], 'Curvature', [1, 1], EdgeColor='k',LineWidth=1, FaceColor=[0,0,0,1]);
rectangle('Position', [d3-radius, H-radius, 2*radius, 2*radius], 'Curvature', [1, 1], EdgeColor='k',LineWidth=1, FaceColor=[0,0,0,1]);

% Add joint parameters to plot
text(d1-width/4, height*1.2, '$$d_1$$', 'Interpreter', 'latex', 'FontSize', 12);
text(d2-width/4, -height*1.2, '$$d_2$$', 'Interpreter', 'latex', 'FontSize', 12);
text(d3-width/4, H+height*1.2, '$$d_3$$', 'Interpreter', 'latex', 'FontSize', 12);

% Draw the End Effector
patch([x1, x2, x3], [y1, y2, y3], 'green', 'edgecolor', 'red') 

% Draw links
plot([d1, x1], [0, y1], 'blue', 'linewidth',2) 
plot([d2, x2], [0, y2], 'blue', 'linewidth',2)
plot([d3, x3], [H, y3], 'blue', 'linewidth',2)

% % ------------------------------------------------ %
% % Singularities add ins
% % MUST add singular direction as input to function
% % ------------------------------------------------ %
% 
% % Draw Singular Direction 
% quiver(x1, y1, singular_direction(1), singular_direction(2),'Linewidth', 1, color='r');
% 
% % Show that links intersect at some point for singular condition
% m1 = (y1-0)/(x1-d1);
% plot([d1-5, x1+5], [0-5*m1, y1+5*m1], 'b--', 'linewidth',1) 
% m2 = (y2-0)/(x2-d2);
% plot([d2-5, x2+5], [0-5*m2, y2+5*m2], 'b--', 'linewidth',1) 
% if x3==d3
%     plot([d3, x3], [H-5, y3+5], 'b--', 'linewidth',1) 
% else
%     m3 = (y3-H)/(x3-d3);
%     plot([d3-5, x3+5], [H-5*m3, y3+5*m3], 'b--', 'linewidth',1) 
% end
% % ------------------------------------------------ %

axis equal
hold off
grid on
axis([-2 5 -1 H+1])
xlabel("x [m]", "Interpreter","latex")
ylabel("y [m]", "Interpreter","latex")
set(gcf, 'Color', 'w'); % set background to white
str = sprintf('x_%.2f_%.2f_%.2f_q_%.2f_%.2f_%.2f', x(1), x(2), x(3), q(1), q(2), q(3));
saveas(gcf, str+".jpg");

end

