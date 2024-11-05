function [H0, C0, G0, H_tilde, C_tilde, G_tilde] = derive_minmax_matrices
    syms l L M m2 m3 theta_1 theta_2 d_3 theta_dot_1 theta_dot_2 d_dot_3 g
    H = [L^2*M + (L^2*m2)/3 + L^2*m3 + M*d_3^2*sin(theta_2)^2 + d_3^2*m3*sin(theta_2)^2 + (l^2*m3*sin(theta_2)^2)/3 - d_3*l*m3*sin(theta_2)^2, -(L*cos(theta_2)*(2*M*d_3 + 2*d_3*m3 - l*m3))/2, -L*sin(theta_2)*(M + m3);
                                                                         -(L*cos(theta_2)*(2*M*d_3 + 2*d_3*m3 - l*m3))/2,    M*d_3^2 + (l^2*m3)/12 + m3*(d_3 - l/2)^2,                    0;
                 -L*sin(theta_2)*(M + m3),                                         0,               M + m3];
    C = [theta_dot_2*(M*d_3^2*cos(theta_2)*sin(theta_2) + d_3^2*m3*cos(theta_2)*sin(theta_2) + (l^2*m3*cos(theta_2)*sin(theta_2))/3 - d_3*l*m3*cos(theta_2)*sin(theta_2)) + d_dot_3*(M*d_3*sin(theta_2)^2 + d_3*m3*sin(theta_2)^2 - (l*m3*sin(theta_2)^2)/2), theta_dot_1*(M*d_3^2*cos(theta_2)*sin(theta_2) + d_3^2*m3*cos(theta_2)*sin(theta_2) + (l^2*m3*cos(theta_2)*sin(theta_2))/3 - d_3*l*m3*cos(theta_2)*sin(theta_2)) - d_dot_3*((L*cos(theta_2)*(2*M + 2*m3))/4 + (L*cos(theta_2)*(M + m3))/2) + (L*theta_dot_2*sin(theta_2)*(2*M*d_3 + 2*d_3*m3 - l*m3))/2, theta_dot_1*(M*d_3*sin(theta_2)^2 + d_3*m3*sin(theta_2)^2 - (l*m3*sin(theta_2)^2)/2) - theta_dot_2*((L*cos(theta_2)*(2*M + 2*m3))/4 + (L*cos(theta_2)*(M + m3))/2);
     - d_dot_3*((L*cos(theta_2)*(2*M + 2*m3))/4 - (L*cos(theta_2)*(M + m3))/2) - theta_dot_1*(M*d_3^2*cos(theta_2)*sin(theta_2) + d_3^2*m3*cos(theta_2)*sin(theta_2) + (l^2*m3*cos(theta_2)*sin(theta_2))/3 - d_3*l*m3*cos(theta_2)*sin(theta_2)),                                                                                                                                                                                                         d_dot_3*(M*d_3 + (m3*(2*d_3 - l))/2),                                 theta_dot_2*(M*d_3 + (m3*(2*d_3 - l))/2) - theta_dot_1*((L*cos(theta_2)*(2*M + 2*m3))/4 - (L*cos(theta_2)*(M + m3))/2);
                                                             theta_dot_2*((L*cos(theta_2)*(2*M + 2*m3))/4 - (L*cos(theta_2)*(M + m3))/2) - theta_dot_1*(M*d_3*sin(theta_2)^2 + d_3*m3*sin(theta_2)^2 - (l*m3*sin(theta_2)^2)/2),                                                                                                                                        theta_dot_1*((L*cos(theta_2)*(2*M + 2*m3))/4 - (L*cos(theta_2)*(M + m3))/2) - theta_dot_2*(M*d_3 + (m3*(2*d_3 - l))/2),                                                                                                                                0];
    
    
    G = [0;     -(g*sin(theta_2)*(2*M*d_3 + 2*d_3*m3 - l*m3))/2;                      g*cos(theta_2)*(M + m3)];
    
    H0 = subs(H, M, 0);
    H_tilde = simplify(H-H0);
    C0 = subs(C, M, 0);
    C_tilde = simplify(C-C0);
    G0 = subs(G, M, 0);
    G_tilde = simplify(G-G0);
end

