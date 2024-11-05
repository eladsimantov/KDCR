function v = v_plan(prof, t)
%V_PLAN Trajectory velocity for different vel profiles
%   prof – decision value for the different velocity profile. 
%   t - time
%   v – the velocity vector of the tool origin in specific time t

pa = [0.4 0 0.8]; pb = [0.25 -0.5 1]; % Hardcoded values
t_i=0; t_f=2; Delta_t=t_f-t_i; % Hardcoded values

% if prof=="Constant"
if prof==1
    v = (pb-pa)/(t_f-t_i);
% elseif prof=="Trapezoidal"
elseif prof==2
    v_max = (6/5)*(pb-pa)/Delta_t;
    if t < Delta_t/6
        v = v_max*(6/Delta_t)*(t);
    elseif t > 5*Delta_t/6
        v = 6*v_max*(1-t/Delta_t);
    else
        v = v_max;
    end
% elseif prof=="Polynomial"
elseif prof==3
    v = (pb-pa) * ((30*t^4)/(Delta_t^5) - (60*t^3)/(Delta_t^4 ) + (30*t^2)/(Delta_t^3));
else
    error("Incorrect Velocity Profile. Use Constant OR Trapezoidal OR Polynomial")
end

if t>=2
    v = [0.0000, 0.0000, 0.0000];
end
end
