function JL = derive_linear_jacobian(A0T, jointVars)
    % A0T: 4x4 homogeneous transformation matrix from the tool to the base
    % jointVars: symbolic variables for the robot's joints 
    % e.g., [theta1, theta2, d3, theta4, theta5]
    
    % Number of joints is the length of the joint variable vector
    nJoints = length(jointVars);
    
    % Initialize the linear Jacobian matrix
    JL = sym(zeros(3, nJoints));
    
    % Loop over each joint variable
    for i = 1:nJoints
        % Compute the partial derivatives of the position part of A0T w.r.t. each joint variable
        JL(1, i) = diff(A0T(1, 4), jointVars(i));
        JL(2, i) = diff(A0T(2, 4), jointVars(i));
        JL(3, i) = diff(A0T(3, 4), jointVars(i));
    end
    
    % Simplify the resulting Jacobian matrix
    JL = simplify(JL);
end
