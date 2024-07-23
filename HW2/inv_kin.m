function [qs] = inv_kin(x)
%INV_KIN Accepts a task vector and returns all possible joint values for
%the Inverse Kinematics solution.
%   x - the task vector of size 3 (x1, y1, theta)
%   theta should be in RAD only.
global r L H

x2 = x(1) + r*cos(x(3));
y2 = x(2) + r*sin(x(3));
x3 = x(1) + r*cos(x(3)+pi/3);
y3 = x(2) + r*sin(x(3)+pi/3);

Delta_1 = L^2 - x(2)^2;
Delta_2 = L^2 - y2^2;
Delta_3 = L^2 - (y3 - H)^2;


if Delta_1<0 
    display(x)
    error("Your task vector results in non exsistant solutions. Please enter a valid task vector")
elseif Delta_2<0
    display(x)
    error("Your task vector results in non exsistant solutions. Please enter a valid task vector")
elseif Delta_3<0
    display(x)
    error("Your task vector results in non exsistant solutions. Please enter a valid task vector")
end

d1_plus = x(1) + (Delta_1)^(1/2);
d1_minus = x(1) - (Delta_1)^(1/2);

d2_plus = x2 + (Delta_2)^(1/2);
d2_minus = x2 - (Delta_2)^(1/2);

d3_plus = x3 + (Delta_3)^(1/2);
d3_minus = x3 - (Delta_3)^(1/2);

qs = [
    [d1_plus d1_plus d1_plus d1_plus d1_minus d1_minus d1_minus d1_minus];
    [d2_plus d2_plus d2_minus d2_minus d2_plus d2_plus d2_minus d2_minus];
    [d3_plus d3_minus d3_plus d3_minus d3_plus d3_minus d3_plus d3_minus]
    ]; % Three joints, Eight solutions.


end


 
 
