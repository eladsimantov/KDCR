function a = a_plan(prof, t)
%A_PLAN Trajectory acceleration for different vel profiles, at any time t.
%   prof – decision value for the different velocity profile. 
%   t - time
%   a – the acceleration of the tool origin in time t

pa = [0.4 0 0.8]; pb = [0.25 -0.5 1]; % Hardcoded values
t_i=0; t_f=2; Delta_t=t_f-t_i; % Hardcoded values for simplicity

% if prof=="Constant"
if prof==1
    a = 0;
% elseif prof=="Trapezoidal"
elseif prof==2
    v_max = (6/5)*(pb-pa)/Delta_t;
    if t < Delta_t/6
        a = v_max*(6/Delta_t);
    elseif t > 5*Delta_t/6
        a = -v_max*(6/Delta_t);
    else
        a = 0;
    end
% elseif prof=="Polynomial"
elseif prof==3
    a = (pb-pa) * ((120*t^3)/(Delta_t^5) - (180*t^2)/(Delta_t^4 ) + (60*t)/(Delta_t^3));
else
    error("Incorrect Velocity Profile. Use Constant OR Trapezoidal OR Polynomial")
end
if t>=2 
    a = [[0.0000, 0.0000, 0.0000]];
end

end
