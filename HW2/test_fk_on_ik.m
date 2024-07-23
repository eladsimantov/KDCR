function [errors] = test_fk_on_ik(x, tol)
%TEST_FK_ON_IK Input a task vector X and check if it is returned in the
%nesting of FK on IK.
%   x - input - task vector
%   tol - input - tolerance for checking equal vectors
%   errors - output - the difference between real task and calculated one.
close all; clc;
display(x)
inv_kin_sols = inv_kin(x); % get all IK solutions for x
display(inv_kin_sols)
N_of_sols = size(inv_kin_sols, 2); % get number of IK solutions 
errors = ones(N_of_sols, 1); % create errors vector per solution
fprintf("Testing IK solution number: \n")
for i=1:N_of_sols
    fprintf(i + " \n")
    fk_sols = forward_kin(inv_kin_sols(1:3,i)); % FK solutions matrix of i-th IK solution 
    err = norm(x - fk_sols(ismembertol(fk_sols, x, tol))); % difference between task x and FK solution
    errors(i) = err;
end
close all;
display(errors)
end

