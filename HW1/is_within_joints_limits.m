function [is_within_limits] = is_within_joints_limits(q)
%IS_JOINTS_WITHIN_LIMITS checks if the joints are within defined limits.
%   q - joint variables as row vector. 
%   Return 1 if joints are within limits.
q1_min = -pi; q1_max = pi;
q2_min = -pi/2; q2_max = pi/2; 
d3_min = 0;
cycle = 2*pi;
% default value
is_within_limits = 0;

% Angles are the most annoying. Here we will deal with possible shifts.
% we will check if the angle, or all possible shifts are within limits. If
% none are withing limits so we know for sure.
% This logic assumes the joint limits are symmetric hence the use of cycle.
% If it is not symmetric simply modify the cases respectiely.

% q(1) logical cases of out of limits
q1_case_1 = and(q(1) <= q1_min, q(1) >= 0 - cycle);
q1_case_2 = and(q(1) <= q1_min + cycle, q(1) >= q1_max);

% q(2) logical cases of out of limits
q2_case_1 = and(q(2) <= q2_min, q(2) >= 0 - cycle);
q2_case_2 = and(q(2) <= q2_min + cycle, q(2) >= q2_max);

if or(q1_case_1, q1_case_2)
    error("q(1) is out of its limit")
elseif or(q2_case_1, q2_case_2)
    error("q(2) is out of its limit")
elseif q(3) <= d3_min
    error("q(3) is out of its limit")
else
    is_within_limits = 1;
end 
end

