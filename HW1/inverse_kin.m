function [q] = inverse_kin(x,elbows)
%INVERSE_KIN returns a single joints solution vector for a given position
%   x - position vector of the tool 
%   elbows – decision values vector of size 1x2 for the different solutions.
%   global parameters - l1 l2 L H r13 r23 r33 r22 r32

% initiallizations
global l1 l2 L H r13 r23 r33
px = x(1);
py = x(2);
pz = x(3);

q = zeros(1,3); % create the row vector 
%q = zeros(1,5); % for the general case of 5 joints 

% in general case l2 is nonzero
px_star = px - l2 * r13;
py_star = py - l2 * r23;
pz_star = pz - l2 * r33;
sqrt_term = (px_star^2 + py_star^2 - L^2);
if sqrt_term < 0
    error("Invalid Point for IK. Try to see that (px_star^2+py_star^2 > L^2) " + ...
        "(px_star = px-l2*r13), (py_star=py-l2*r23)")
end

% theta1 and theta2 solutions by elbows
if elbows==[1 1]
    % Solve for theta1 elbow up and theta2 upper
    q(1) = atan2(py_star, px_star) + atan2(sqrt(sqrt_term), L); % theta1
    s2 = (px_star - L*cos(q(1))) / sin(q(1));
    c2 = pz_star - H;
    q(2) = atan2(s2, c2); % theta2
elseif elbows==[-1 2]
    % Solve for theta1 elbow down and theta2 upper
    q(1) = atan2(py_star, px_star) + atan2(-sqrt(sqrt_term), L); % theta1
    s2 = (px_star - L*cos(q(1))) / sin(q(1));
    c2 = pz_star - H;
    q(2) = atan2(s2, c2); % theta2
elseif elbows==[1 -1]
    % Solve for theta1 elbow up and theta2 lower
    q(1) = atan2(py_star, px_star) + atan2(sqrt(sqrt_term), L); % theta1
    s2 = (px_star - L*cos(q(1))) / sin(q(1));
    c2 = pz_star - H;
    q(2) = atan2(-s2, -c2); % theta2 lower
elseif elbows==[-1 -1]
    % Solve for theta1 elbow down and theta2 lower
    q(1) = atan2(py_star, px_star) + atan2(-sqrt(sqrt_term), L); % theta1
    s2 = (px_star - L*cos(q(1))) / sin(q(1));
    c2 = pz_star - H;
    q(2) = atan2(-s2, -c2); % theta2 lower 
else
    error('Elbow indexing incorrect. Use [1 1], [1 -1] [-1 1] and [-1 -1] for all four cases.');
end

% % Solve for theta2
% s2 = (px_star - L*cos(q(1))) / sin(q(1));
% c2 = pz_star - H;
% q(2) = atan2(s2, c2);
% % q(2) = atan2(-s2, -c2) another case for negative d3 but not physical

% Solve for d3
q(3) = (pz_star - H) / cos(q(2)) - l1;

% Find theta4 generally for future use
% s2=sin(q(2)); c2=cos(q(2)); % Just for notations
% s4 = (-s2*r22 + cos(q(1))*c2*r32)/(sin(q(1))*s2);
% c4 = r32 / s2;
% q(4) = atan2(s4, c4);
% % Find theta5 generally for future use
% s5 = (-c2*r13 + s2*s4*r33)/(-s2^2*s4^2-c2^2);
% c5 = (-s2*s4*r13 - c2*r33)/(-s2^2*s4^2-c2^2);
% q(5) = atan2(s5, c5); 
end


