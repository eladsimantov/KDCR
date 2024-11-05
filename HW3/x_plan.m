function x = x_plan(prof,t)
% This function receives a decision value for the velocity profile and
% a time 𝑡 and returns the tool position.
t1 = 1/3;
t2 = 5/3;
x0 = 0.4;
y0 = 0;
z0 = 0.8;
ax = -0.27;
ay = -0.9;
az = 0.36;
x1 = x0+ax*t1^2*0.5;
y1 = y0+ay*t1^2*0.5;
z1 = z0+az*t1^2*0.5;
x2 = x1+ax*t1*(t2-t1);
y2 = y1+ay*t1*(t2-t1);
z2 = z1+az*t1*(t2-t1);
switch prof
    case 1
        x = [-0.075*t+x0,-0.25*t+y0,0.1*t+z0];
    case 2
        if t<t1
            x = [x0+ax*t^2*0.5 y0+ay*t^2*0.5 z0+az*t^2*0.5];
        elseif t<t2
            
            x = [x1+ax*t1*(t-t1) y1+ay*t1*(t-t1) z1+az*t1*(t-t1)];
        else
            x = [x2+ax*t1*(t-t2)-ax*(t-t2)^2*0.5 y2+ay*t1*(t-t2)-ay*(t-t2)^2*0.5...
                z2+az*t1*(t-t2)-az*(t-t2)^2*0.5];
        end
    case 3
        x = [0.4-0.1875*t^3+0.1406*t^4-0.0281*t^5,...
            -0.625*t^3+0.4688*t^4-0.0938*t^5,...
            0.8+0.25*t^3-0.1875*t^4+0.0375*t^5];
end
end