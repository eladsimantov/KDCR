function  tau=tau_cont(q,qdot,t,law) 
%TAU_CONT choose the control law "ID" or "G" or "PID" or "MINMAX"


prof = 3; % HARDCODED... Polynomial profile.
qdes        =   q_plan(prof,t);                 % Desired q
qdes_dot    =   q_dot_plan(prof,t);             % Desired q_dot
qdes_dot2   =   q_dot2_plan(prof,t);            % Desired q_dot2

switch law
    case "ID"
        % The ID control law includes the proportional Kp and derivative Kd
        % parameters along with precisely measured q, qdot values.
        M = 0.5;
        % Bessel prototype of 2nd order for minimal oscillation and fast response
        w0 = 4; % for sufficient settling time, and small for minimal control effort
        Kp = w0^2; 
        Kd = 2*0.866*w0; 
        
        [H_mat,C_mat,G] = dynamics_mat(q,qdot, M);
        tau = C_mat*qdot' + G + H_mat*(qdes_dot2' - Kp*(q' - qdes') - Kd*(qdot' - qdes_dot'));
    case "G"
        % The G control law only has vector G known in the dynamics mat
        M = 0.5;
        % Set Symmetric and Positive definite for set point local stability 
        Kp = [1000 0 0; 0 1200 0; 0 0 1100]; 
        Kd = [30 0 0; 0 40 0; 0 0 50]; 
        [~,~,G] = dynamics_mat(q, qdot, M);
        tau = G - Kp*(q' - qdes') - Kd*(qdot' - qdes_dot'); % Make sure all vector sizes are correct

    case "PID"
        % 0<M<1  Mass is randomly chosen in main according to rng(1).    
        global KpPID KiPID KdPID X_numPID
        % This is not ideal but the project criteria forced us to use global because of the requirement for "tau_cont" function format. 
        % some may suggest "persistent" but that is not a good option because the value for eI must be shared between "tau_cont" and "state_eq_PID" functions, 
        % and since the ode45 does runga cutta we cannot simply calculate eI in the same way as done by ode45 !!!
        e       = q - qdes;
        edot    = qdot - qdes_dot;
        idx     = 1 + round(t/0.001);
        eI      = X_numPID(idx, 1:3);
        tau     = -KpPID*e' - KdPID*edot' - KiPID*eI';


    case "MINMAX"
        % 0<M<1  Mass is randomly chosen in main according to rng(2).
        % define MINMAX parameters
        K       = [50 0 0; 0 30 0; 0 0 50];
        P       = [30 0 0; 0 30 0 ; 0 0 30];
        beta    = 1.4;
        delta   = 0.005;

        e       = q' - qdes';           % As defined in equation (4.2) in WORD
        edot    = qdot' - qdes_dot';    % As defined in equation (4.2) in WORD
        s       = edot + K*e;           % As defined in equation (4.3) in WORD
        qdot_r  = qdot' - s;            % As defined in equation (4.4) in WORD
        qddot_r = qdes_dot2' - K*edot;  % As defined in equation (4.5) in WORD
        
        [H0,C0,G0] = dynamics_mat(q,qdot,0);
        [H1,C1,G1] = dynamics_mat(q,qdot,1);
        eta0 = -H0*qddot_r-C0*qdot_r-G0+P*e;    % As defined in equation (4.12) in WORD

        % Finding ρ̃  as defined in equation (4.16) in WORD 
        rho_tilde = beta * (norm((H1-H0) * qddot_r) + norm((C1-C0) * qdot_r) + norm(G1-G0));
        
        % Finding ρ_u as defined in equation (4.17) in WORD
        rho_u = max(0, 1/(norm(s)+delta)*(s'*eta0-e'*K'*P*e) + rho_tilde); 
        tau = -rho_u*s/(norm(s)+delta); % As defined in equation (4.18) in WORD

end

