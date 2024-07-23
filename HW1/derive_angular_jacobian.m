function JA = derive_angular_jacobian(A0T_as_func_of_t, joint_symbols_as_func_of_t, t_symbol)
    % A0T: 4x4 homogeneous transformation matrix from the tool to the base
    % as a function of t.
    % joint_symbols_as_func_of_t: symbolic variables for the robot's joints
    % as a function of t. e.g, [theta_1(t) theta_2(t) d_3(t) theta_4(t)].
    
    % Number of joints is the length of the joint variable vector
    nJoints = length(joint_symbols_as_func_of_t);
    
    % Initialize the angular Jacobian matrix
    JA = sym(zeros(3, nJoints));
    
    % Calc the rotation matrix, its derivative and then Omega
    R0T = A0T_as_func_of_t(1:3,1:3);
    R0T_dot = diff(R0T, t_symbol);
    Omega = R0T_dot * transpose(R0T);
    Omega_Simplified = simplify(Omega);
    q_dot = diff(joint_symbols_as_func_of_t, t_symbol);

    % Extract the omega vector from the Omega Matrix
    omega_vec = [Omega_Simplified(3,2); Omega_Simplified(1,3); Omega_Simplified(2,1)];

    % Loop over each joint variable
    for i = 1:nJoints
        % Compute the partial derivatives
        JA(1, i) = diff(omega_vec(1), q_dot(i));
        JA(2, i) = diff(omega_vec(2), q_dot(i));
        JA(3, i) = diff(omega_vec(3), q_dot(i));
    end
    
    % Simplify the resulting Jacobian matrix
    JA = simplify(JA, "Steps", 100);
end
