function [x] = forward_kin(q)
%FORWARD_KIN
%   q – row vector of the joints parameters
%   x – position row vector of the tool
global l1 l2 L H

% check limitations on joint variables.
is_within_joints_limits(q);
x = [L*cos(q(1)) + q(3)*sin(q(1))*sin(q(2)) + l1*sin(q(1))*sin(q(2)) + l2*sin(q(1))*sin(q(2)); 
    L*sin(q(1)) - q(3)*cos(q(1))*sin(q(2)) - l1*cos(q(1))*sin(q(2)) - l2*cos(q(1))*sin(q(2));
    H + q(3)*cos(q(2)) + l1*cos(q(2)) + l2*cos(q(2))]';
end

