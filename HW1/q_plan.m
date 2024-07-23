function [Q] = q_plan(prof, t, elbows)
%Q_PLAN Takes a prof and a time vector t, and 
%returns the joint parameters Matrix at all
%times t by the use of IK and x_plan function.
%Q has 3 columns for 3 joints and one row per value of t.
%   t - the time vector.
%   elbows - solution decision variables vector used in the IK solver.
%   Q - returned value - the joints parameters matrix in all times t.
Q = inverse_kins(x_plan(prof, t), elbows);
end

