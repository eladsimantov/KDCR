function [V] = v_plan(prof, t)
%V_PLAN Trajectory velocity for different vel profiles
%   prof – decision value for the different velocity profile. 
%   t - times vector where i=length(t)
%   v – the velocity vector of the tool origin in specific time t(i).
%   V - returned value - the velocity Matrix of the tool in all times t.

V = zeros(length(t), 3); % init
pa = [0.4 0 0.8]; pb = [0.25 -0.5 1]; % Hardcoded values
t_i=0; t_f=2; Delta_t=t_f-t_i; % Hardcoded values
for i=1:length(t)

    if prof=="Constant"
        v = (pb-pa)/(t_f-t_i);
    elseif prof=="Trapezoidal"
        v_max = (6/5)*(pb-pa)/Delta_t;
        if t(i) < Delta_t/6
            v = v_max*(6/Delta_t)*(t(i));
        elseif t(i) > 5*Delta_t/6
            v = 6*v_max*(1-t(i)/Delta_t);
        else
            v = v_max;
        end
    elseif prof=="Polynomial"
        v = (pb-pa) * ((30*t(i)^4)/(Delta_t^5) - (60*t(i)^3)/(Delta_t^4 ) + (30*t(i)^2)/(Delta_t^3));
    else
        error("Incorrect Velocity Profile. Use Constant OR Trapezoidal OR Polynomial")
    end
    V(i, :) = v; % insert velocity at time t(i)
end 

end
