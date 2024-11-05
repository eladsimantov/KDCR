function [H_mat,C_mat,G] = dynamics_mat(q,q_dot, M)
    global H L l
    d       = 0.015;
    rho     = 7800;
    m1      = pi*(d/2)^2*H*rho;
    m2      = pi*(d/2)^2*L*rho;
    m3      = pi*(d/2)^2*l*rho;
    % M       = 0.5;
    % M       = 0; % Run simulation with no load mass for Q3
    g       = 9.81;
    th1     = q(1);
    th2     = q(2);
    d3      = q(3);
    th1_d   = q_dot(1);
    th2_d   = q_dot(2);
    d3_d    = q_dot(3);
    H_mat   = [L^2*M + (L^2*m2)/3 + L^2*m3 + M*d3^2*sin(th2)^2 + d3^2*m3*sin(th2)^2 + (l^2*m3*sin(th2)^2)/3 - d3*l*m3*sin(th2)^2, -(L*cos(th2)*(2*M*d3 + 2*d3*m3 - l*m3))/2, -L*sin(th2)*(M + m3);
                                                                        -(L*cos(th2)*(2*M*d3 + 2*d3*m3 - l*m3))/2,    M*d3^2 + (l^2*m3)/12 + m3*(d3 - l/2)^2,                    0;
                                                                                             -L*sin(th2)*(M + m3),                                         0,               M + m3];
 
    C_mat   = [th2_d*(M*d3^2*cos(th2)*sin(th2) + d3^2*m3*cos(th2)*sin(th2) + (l^2*m3*cos(th2)*sin(th2))/3 - d3*l*m3*cos(th2)*sin(th2)) + d3_d*(M*d3*sin(th2)^2 + d3*m3*sin(th2)^2 - (l*m3*sin(th2)^2)/2), th1_d*(M*d3^2*cos(th2)*sin(th2) + d3^2*m3*cos(th2)*sin(th2) + (l^2*m3*cos(th2)*sin(th2))/3 - d3*l*m3*cos(th2)*sin(th2)) - d3_d*((L*cos(th2)*(2*M + 2*m3))/4 + (L*cos(th2)*(M + m3))/2) + (L*th2_d*sin(th2)*(2*M*d3 + 2*d3*m3 - l*m3))/2, th1_d*(M*d3*sin(th2)^2 + d3*m3*sin(th2)^2 - (l*m3*sin(th2)^2)/2) - th2_d*((L*cos(th2)*(2*M + 2*m3))/4 + (L*cos(th2)*(M + m3))/2);
                - d3_d*((L*cos(th2)*(2*M + 2*m3))/4 - (L*cos(th2)*(M + m3))/2) - th1_d*(M*d3^2*cos(th2)*sin(th2) + d3^2*m3*cos(th2)*sin(th2) + (l^2*m3*cos(th2)*sin(th2))/3 - d3*l*m3*cos(th2)*sin(th2)),                                                                                                                                                                                                         d3_d*(M*d3 + (m3*(2*d3 - l))/2),                                 th2_d*(M*d3 + (m3*(2*d3 - l))/2) - th1_d*((L*cos(th2)*(2*M + 2*m3))/4 - (L*cos(th2)*(M + m3))/2);
                th2_d*((L*cos(th2)*(2*M + 2*m3))/4 - (L*cos(th2)*(M + m3))/2) - th1_d*(M*d3*sin(th2)^2 + d3*m3*sin(th2)^2 - (l*m3*sin(th2)^2)/2),                                                                                                                                        th1_d*((L*cos(th2)*(2*M + 2*m3))/4 - (L*cos(th2)*(M + m3))/2) - th2_d*(M*d3 + (m3*(2*d3 - l))/2),                                                                                                                                0];
 
    G       = [0.0;-(g*sin(th2)*(2*M*d3 + 2*d3*m3 - l*m3))/2;g*cos(th2)*(M + m3)];
end

