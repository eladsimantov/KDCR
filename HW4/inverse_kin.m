function q=inverse_kin(x,elbows,R)
%  This function receives the position vector of the tool,
% the decision values vector for the different solutions, and the rotation of the tool.
% The variable 𝑞 represents the joint values.
% If 𝑅 is not a matrix, the simulation is run.
global H L;
l1 = 0;
l2 = 0;
q = zeros(size(x,1),3);
for i=1:size(q,1)
    px = x(i,1);
    py = x(i,2);
    pz = x(i,3);
    alpha = sqrt(px^2+py^2-L^2);
    d3 = elbows(1)*sqrt(alpha^2+(pz-H)^2)-l1;
    theta2 = atan2(elbows(2)*alpha/(d3+l1),(pz-H)/(d3+l1));
    theta1 = atan2((alpha*px+L*py),(L*px-alpha*py));
    if size(R) == 3
        ryx = R(1,2);
        ryy = R(2,2);
        ryz = R(3,2);
        if sin(theta1) == 0
            theta4 = atan2(-cos(theta1)*cos(theta2)*ryy,cos(theta1)*ryx);
        elseif sin(theta2) == 0
            theta4 = atan2(-cos(theta1)*cos(theta2)*ryx-...
                cos(theta2)*sin(theta1)*ryy,...
                -sin(theta1)*ryx+cos(theta1)*ryy);
        else
            theta4 = atan2((cos(theta1)*cos(theta2)*ryz-sin(theta2)*ryy)/...
                (sin(theta1)*sin(theta2)),ryz/sin(theta2));
        end
        theta5 = atan2(-cos(theta2)*rxz+sin(theta2)*sin(theta4)*rzz,...
            -sin(theta2)*sin(theta4)*txz-cos(theta2)*rzz);
    else
        theta4 = 0;
        theta5 = 0;
    end
    q(i,:) = [theta1,theta2,d3];
end
end