function a = a_plan(prof,t)
% This function receives a decision value for the velocity profile
% and a time 𝑡, and returns the tool acceleration.
t1 = 1/3;
t2 = 5/3;
ax = -0.27;
ay = -0.9;
az = 0.36;

switch prof
    case 1
        a = [0,0,0];
    case 2
        if t<t1
            a = [ax ay az];
        elseif t<t2
            
            a = [0 0 0];
        else
            a = [-ax -ay...
                -az];
        end
    case 3
        a = [-0.1875*6*t+0.1406*12*t^2-0.0281*20*t^3,...
            -0.625*6*t+0.4688*12*t^2-0.0938*20*t^3,...
            0.25*6*t-0.1875*12*t^2+0.0375*20*t^3];
end
end