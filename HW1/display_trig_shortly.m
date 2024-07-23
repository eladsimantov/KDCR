function displayed_object = display_trig_shortly(symbolic_object, symbolic_joint_angles, symbolic_replaced_values)
%DISPLAY_TRIG_SHORTLY takes a symbolic object and replaces all the cosine and sine functions with c and s for short. 
    % symbolic_object: e.g., a Matrix A0T
    % joint_params: the joint angles e.g., [theta2 theta5]
    % symbolic_replaced_values: cosine in top row, sin in bottom. e.g., [c2 c5 ;s2 s5]
    % Number of joints is the length of the joint variable vector
    nAngles = length(symbolic_joint_angles);
    
    % Initialize the sine and cosine vector
    old_vec = sym(zeros(2, nAngles));

    % Loop over each joint angle to replace its cosine and sine
    for i = 1:nAngles
        old_vec(1, i) = cos(symbolic_joint_angles(i)); % enter the cosines 
        old_vec(2, i) = sin(symbolic_joint_angles(i)); % enter the sines
    end 
    % 
    displayed_object = subs(symbolic_object, old_vec, symbolic_replaced_values);
end

