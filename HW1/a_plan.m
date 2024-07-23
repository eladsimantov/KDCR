function [A] = a_plan(prof, t)
%A_PLAN Trajectory acceleration for different vel profiles, at all times t.
%   prof – decision value for the different velocity profile. 
%   t - times vector where i=length(t)
%   a – the acceleration of the tool origin in time t(i)
%   A - returned value - the acceleration Matrix of the tool in all times t

A = zeros(length(t), 3); % init
pa = [0.4 0 0.8]; pb = [0.25 -0.5 1]; % Hardcoded values
t_i=0; t_f=2; Delta_t=t_f-t_i; % Hardcoded values for simplicity

for i = 1:length(t)
    if prof=="Constant"
        a = 0;
    elseif prof=="Trapezoidal"
        v_max = (6/5)*(pb-pa)/Delta_t;
        if t(i) < Delta_t/6
            a = v_max*(6/Delta_t);
        elseif t(i) > 5*Delta_t/6
            a = -v_max*(6/Delta_t);
        else
            a = 0;
        end
    elseif prof=="Polynomial"
        a = (pb-pa) * ((120*t(i)^3)/(Delta_t^5) - (180*t(i)^2)/(Delta_t^4 ) + (60*t(i))/(Delta_t^3));
    else
        error("Incorrect Velocity Profile. Use Constant OR Trapezoidal OR Polynomial")
    end
    A(i, :) = a; % insert accel at time t(i) into A matrix.

end

end
