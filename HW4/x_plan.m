function x = x_plan(prof,t)
%X_PLAN Trajectory locator for different vel profiles that returns a
%position Matrix X for any time t.
%   prof – decision value for the different velocity profile. 
%   t - time
%   x – the position row vector of the tool origin in speicific time t.

pa = [0.4 0 0.8]; pb = [0.25 -0.5 1]; % Hardcoded values
t_i=0; t_f=2; Delta_t=t_f-t_i; % Hardcoded values for simplicity

% if prof=="Constant"
if prof==1
    x = pa + (pb-pa)*(t)/Delta_t;
% elseif prof=="Trapezoidal"
elseif prof==2
    if t < Delta_t/6
        x = pa + (pb-pa) * (36/10) * ((t)/Delta_t)^2;
    elseif t > 5*Delta_t/6
        x = pa + (pb-pa)*(9/10) + (pb-pa)*(36/5)*(t/Delta_t)*(1 - t/(2*Delta_t)) - (pb-pa)*(7/2);
    else
        x = pa + (pb-pa)/10 + (6/5)*(pb-pa)*((t/Delta_t) - 1/6);
    end
% elseif prof=="Polynomial"
elseif prof==3
    x = (pb-pa) * ((6*t^5)/(Delta_t^5) - (15*t^4)/(Delta_t^4 ) + (10*t^3)/(Delta_t^3 )) + pa;
else
    error("Incorrect Velocity Profile. Use Constant OR Trapezoidal OR Polynomial")
end
 
if t>=2
    x = pb;
end

end

