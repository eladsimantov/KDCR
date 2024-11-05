function q=q_plan(prof,t) 
% This function receives a decision value for the velocity profile and the time t 
% and returns the joints valuse
x = x_plan(prof,t);
global elbow1 elbow2 elbow3
R = 1;
elbows = [elbow1 elbow2 elbow3];
q = inverse_kin(x,elbows,R);
if q(3) < 0
    elbow1 = -1*elbow1;
    q = inverse_kin(x,elbows,R);
end
if abs(q(1))>pi
    elbow3 = -1*elbow3;
    q = inverse_kin(x,elbows,R);
end
if abs(q(2))>pi/2
    elbow2 = -1*elbow2;
    q = inverse_kin(x,elbows,R);
end

end