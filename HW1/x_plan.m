function [X] = x_plan(prof,t)
%X_PLAN Trajectory locator for different vel profiles that returns a
%position Matrix X for all times t.
%   prof – decision value for the different velocity profile. 
%   t - times vector where i=length(t)
%   x – the position row vector of the tool origin in speicific time t.
%   X - the position matrix of the tool origin in all times t.

X = zeros(length(t), 3); % init
pa = [0.4 0 0.8]; pb = [0.25 -0.5 1]; % Hardcoded values
t_i=0; t_f=2; Delta_t=t_f-t_i; % Hardcoded values for simplicity

for i=1:length(t)
    if prof=="Constant"
        x = pa + (pb-pa)*(t(i))/Delta_t;
    elseif prof=="Trapezoidal"
        if t(i) < Delta_t/6
            x = pa + (pb-pa) * (36/10) * ((t(i))/Delta_t)^2;
        elseif t(i) > 5*Delta_t/6
            x = pa + (pb-pa)*(9/10) + (pb-pa)*(36/5)*(t(i)/Delta_t)*(1 - t(i)/(2*Delta_t)) - (pb-pa)*(7/2);
        else
            x = pa + (pb-pa)/10 + (6/5)*(pb-pa)*((t(i)/Delta_t) - 1/6);
        end
    elseif prof=="Polynomial"
        x = (pb-pa) * ((6*t(i)^5)/(Delta_t^5) - (15*t(i)^4)/(Delta_t^4 ) + (10*t(i)^3)/(Delta_t^3 )) + pa;
    else
        error("Incorrect Velocity Profile. Use Constant OR Trapezoidal OR Polynomial")
    end
    X(i, :) = x;
end 

end

