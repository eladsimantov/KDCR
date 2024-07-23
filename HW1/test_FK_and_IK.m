function [result] = test_FK_and_IK(test, vec)
%TEST_IK_AND_FK will test the IK and FK. enter a vector and a test type and
%get a T/F result on your test.
%   test: Enter 0 OR 1 for IK on FK or FK on IK tests.
%   vec: Enter q for test=0 or x for test=1
if test==1
    q = vec;
    calculated_q_up = inverse_kin(forward_kin(q), 1);
    calculated_q_down = inverse_kin(forward_kin(q), -1);
    result = min(norm(calculated_q_down - q), norm(calculated_q_up - q));
else
    x = vec;
    calculated_x_elbow_up = forward_kin(inverse_kin(x, 1));
    calculated_x_elbow_down = forward_kin(inverse_kin(x, -1));
    result = min(norm(calculated_x_elbow_up - x), norm(calculated_x_elbow_down - x));
end
end


