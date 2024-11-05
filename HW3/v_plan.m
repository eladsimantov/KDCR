function v = v_plan(prof,t)
% This function receives a decision value for the velocity profile and
% a time t, and returns the tool velocity.
t1 = 1/3;
t2 = 5/3;
ax = -0.27;
ay = -0.9;
az = 0.36;

switch prof
    case 1
        v = [-0.075,-0.25,0.1];
    case 2
        if t<t1
            v = [ax*t ay*t az*t];
        elseif t<t2
            
            v = [ax*t1 ay*t1 az*t1];
        else
            v = [ax*t1-ax*(t-t2) ay*t1-ay*(t-t2)...
                az*t1-az*(t-t2)];
        end
    case 3
        v = [-0.1875*3*t^2+0.1406*4*t^3-0.0281*5*t^4,...
            -0.625*3*t^2+0.4688*4*t^3-0.0938*5*t^4,...
            0.25*3*t^2-0.1875*4*t^3+0.0375*5*t^4];
end
end