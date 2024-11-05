%% Main Project 3 code
%% Symbolic section
clear all
global H L l
syms H L g l real
syms q1 q2 q3 r1 r2 r3 real
syms c1 s1 c2 s2
syms q_dot_1 q_dot_2 q_dot_3 real
syms I1 I2 I3 I_t real
syms m1 m2 m3 M real
trig_arr = [c1 s1 c2 s2];
cossin_arr = [cos(q1) sin(q1) cos(q2) sin(q2)];
q = [q1,q2,q3];
q_dot = [q_dot_1,q_dot_2,q_dot_3];
x = sym('x',[size(q,1) 3]);
% Calulating x:
for i=1:size(x,1)
    % The transformation matrixs:
    theta_1 = q(i,1);
    theta_2 = q(i,2);
    d_3 = q(i,3);
    A01 = [...
        cos(theta_1) -sin(theta_1) 0 0;...
        sin(theta_1) cos(theta_1) 0 0;...
        0 0 1 H;...
        0 0 0 1];
    A12 = [...
        1 0 0 L;...
        0 cos(theta_2) -sin(theta_2) 0;...
        0 sin(theta_2) cos(theta_2) 0;...
        0 0 0 1];
    A23 = [...
        1 0 0 0;...
        0 1 0 0;...
        0 0 1 d_3;...
        0 0 0 1];
    
end
A02 = A01*A12;
A0t = A02*A23;

% Rotation matrixes
R01 = A01(1:3,1:3);
R02 = A02(1:3,1:3);
R0t = A0t(1:3,1:3);

% Origen location world frame
Ori_1 = A01(1:3,4);
Ori_2 = A02(1:3,4);
Ori_t = A0t(1:3,4);

% World frame vectors
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Center of mass in linck frame
r1 = H/2;
r2 = L/2;
r3 = l/2;
rc1 = r1*e3;
rc2 = r2*e1;
rc3 = (q3-r3)*e3;

% Center of mass in world frame
Rc1 = rc1;
Rc2 = R01*rc2+Ori_1;
Rc3 = R02*rc3+Ori_2;
Rct = Ori_t;

% linear velocite
Rc1_dot = jacobian(Rc1,q)*q_dot.';
Rc2_dot =  jacobian(Rc2,q)*q_dot.';
Rc3_dot =  jacobian(Rc3,q)*q_dot.';
Rct_dot =  jacobian(Rct,q)*q_dot.';

% moment of inertia in link frame
I1 = (m1*H^2)/12;
I2 = (m2*L^2)/12;
I3 = (m3*l^2)/12;
I_t = M*q3^2;
I_tilde_1 = I1*[1,0,0;
    0,1,0;
    0,0,0];
I_tilde_2 = I2*[0,0,0;
    0,1,0;
    0,0,1];
I_tilde_3 = I3*[1,0,0;
    0,1,0;
    0,0,0];

I_tilde_t = I_t*[1,0,0;
    0,1,0;
    0,0,0];
% moment of inertia in world frame
Inertia1 = I_tilde_1;
Inertia2 = R01*I_tilde_2*R01';
Inertia3 = R02*I_tilde_3*R02';
Inertiat = R02*I_tilde_t*R02';
% Angular velocity of the links
omega1 = q_dot_1*e3;
omega2 = omega1 + q_dot_2*(R01*e1);
omega3 = omega2;


%Q1.a
% Potrntioal Enargy
V1 = m1*g*dot(Rc1,e3);
V2 = m2*g*dot(Rc2,e3);
V3 = m3*g*dot(Rc3,e3);
Vt = M*g*dot(Rct,e3);
V_tot = V1 + V2 + V3 + Vt;

% Kinetick energy
T1 = 0.5*m1*Rc1_dot'*Rc1_dot+0.5*omega1'*Inertia1*omega1;
T2 = 0.5*m2*Rc2_dot'*Rc2_dot+0.5*omega2'*Inertia2*omega2;
T3 = 0.5*m3*Rc3_dot'*Rc3_dot+0.5*omega3'*Inertia3*omega3;
Tt = 0.5*M*Rct_dot'*Rct_dot+0.5*omega3'*Inertiat*omega3;
T_tot = T1 + T2 + T3 + Tt;

syms tau1 tau2 f3 real
deltaW_jonts = [tau1;tau2;f3];
% end effector
syms F1 F2 F3 real
Fe = [F1;F2;F3];
W_Fe = dot(Fe,Ori_t);
deltaW_Fe = jacobian(W_Fe,q)';
% Calculating H and G:
H_mat = jacobian(jacobian(T_tot,q_dot)',q_dot);
G_mat = jacobian(V_tot,q)';
%% Non-con forces %%
% Joints
syms tau1 tau2 f3 real
deltaW_jonts = [tau1;tau2;f3];
% end effector
syms F1 F2 F3 real
Fe = [F1;F2;F3];
W_Fe = dot(Fe,Ori_t);
deltaW_Fe = jacobian(W_Fe,q)';

% Calculating H and G:
H_mat = jacobian(jacobian(T_tot,q_dot)',q_dot);
G_mat = jacobian(V_tot,q)';

%Q1.b
% Calculating links jacobians
J = jacobian_mat_simbolic(q); % Jacobian in the world frame
J_T = [R0t'*J(1:3,:);R0t'*J(4:6,:)];
J1 = [subs(J(:,1),[H,L,l,q2,q3],[H/2,0,0,0,0]),zeros(6,2)];
J2 = [subs(J(:,1:2),[L,l,q3],[L/2,0,0]),zeros(6,1)];
J3 = subs(J,q3,q3-l/2);
% J1 = [subs(J_T(:,1),[H,L,l,q2,q3],[H/2,0,0,0,0]),zeros(6,2)];
% J2 = [subs(J_T(:,1:2),[L,l,q3],[L/2,0,0]),zeros(6,1)];
% J3 = subs(J_T,q3,q3-l/2);
G_b = g*(m1*J1(1:3,:)'+m2*J2(1:3,:)'+(m3)*J3(1:3,:)'+M*J(1:3,:)')*e3;
H_b = m1*J1(1:3,:)'*J1(1:3,:)+m2*J2(1:3,:)'*J2(1:3,:)+(m3)*J3(1:3,:)'*J3(1:3,:) + M*J(1:3,:)'*J(1:3,:)+ J1(4:6,:)'*Inertia1*J1(4:6,:) + J2(4:6,:)'*Inertia2*J2(4:6,:)+J3(4:6,:)'*Inertia3*J3(4:6,:)+J(4:6,:)'*Inertiat*J(4:6,:);
% Calculating C
C = sym(zeros(3,3));
c = sym(zeros(1,3));
for i = 1:3
    for j = 1:3
        for k = 1:3
            c(k) = 0.5*(diff(H_b(i,j),q(k))+diff(H_b(i,k),q(j))-diff(H_b(j,k),q(i)));
        end
        C(i,j) = c*q_dot';
    end
end
%% Numeric calculaition
%% Q2
clear all
global H L l elbow1 elbow2 elbow3
H = 0.2;
L = 0.1;
l = 0.6;
elbow1 = 1;
elbow2 = 1;
elbow3 = 1;

% plotting for prof = 3
prof = 3;
t = 0:0.001:2;
tau = [];
for i = t
    tau = [tau ,tau_plan(prof,i)];
end
figure();
h = tiledlayout(3,1);
nexttile
plot(t,tau(1,:),'DisplayName','\tau_1')
legend
grid on
ylabel('\tau_1 [Nm]')
set(gca,'FontSize',14);
nexttile
plot(t,tau(2,:),'DisplayName','\tau_2')
legend
grid on
ylabel('\tau_2 [Nm]')
set(gca,'FontSize',14);
nexttile
plot(t,tau(3,:),'DisplayName','F{_3}')
ylabel('F[N]')
xlabel('Time[s]')
legend
grid on
title(h,'Torques and forces at the joints for prof 3')
set(gca,'FontSize',14);
% % prof = 2;
% % t = 0:0.001:2;
% % tau = [];
% % for i = t
% %     tau = [tau ,tau_plan(prof,i)];
% % end
% % figure();
% % h = tiledlayout(3,1);
% % nexttile
% % plot(t,tau(1,:),'DisplayName','\tau_1')
% % legend
% % grid on
% % ylabel('\tau_1 [Nm]')
% % set(gca,'FontSize',14);
% % nexttile
% % plot(t,tau(2,:),'DisplayName','\tau_2')
% % legend
% % grid on
% % ylabel('\tau_2 [Nm]')
% % set(gca,'FontSize',14);
% % nexttile
% % plot(t,tau(3,:),'DisplayName','F{_3}')
% % ylabel('F[N]')
% % xlabel('Time[s]')
% % legend
% % grid on
% % title(h,'Torques and forces at the joints for prof 2')
% % set(gca,'FontSize',14);
% % prof = 1;
% % t = 0:0.001:2;
% % tau = [];
% % for i = t
% %     tau = [tau ,tau_plan(prof,i)];
% % end
% % figure();
% % h = tiledlayout(3,1);
% % nexttile
% % plot(t,tau(1,:),'DisplayName','\tau_1')
% % legend
% % grid on
% % ylabel('\tau_1 [Nm]')
% % set(gca,'FontSize',14);
% % nexttile
% % plot(t,tau(2,:),'DisplayName','\tau_2')
% % legend
% % grid on
% % ylabel('\tau_2 [Nm]')
% % set(gca,'FontSize',14);
% % nexttile
% % plot(t,tau(3,:),'DisplayName','F{_3}')
% % ylabel('F[N]')
% % xlabel('Time[s]')
% % legend
% % grid on
% % set(gca,'FontSize',14);
% % title(h,'Torques and forces at the joints for prof 1')
% % set(gca,'FontSize',14);
%% Q3 simbolic
global H L l
syms H L g l real
syms q1 q2 q3 r1 r2 r3 real
syms c1 s1 c2 s2
syms q_dot_1 q_dot_2 q_dot_3 real
syms q_ddot_1 q_ddot_2 q_ddot_3 real
syms I1 I2 I3 I_t real
syms m1 m2 m3 M real
trig_arr = [c1 s1 c2 s2];
cossin_arr = [cos(q1) sin(q1) cos(q2) sin(q2)];
q = [q1,q2,q3];
q_dot = [q_dot_1,q_dot_2,q_dot_3];
x = sym('x',[size(q,1) 3]);
% Calulating x:
for i=1:size(x,1)
    % The transformation matrixs:
    theta_1 = q(i,1);
    theta_2 = q(i,2);
    d_3 = q(i,3);
    A01 = [...
        cos(theta_1) -sin(theta_1) 0 0;...
        sin(theta_1) cos(theta_1) 0 0;...
        0 0 1 H;...
        0 0 0 1];
    A12 = [...
        1 0 0 L;...
        0 cos(theta_2) -sin(theta_2) 0;...
        0 sin(theta_2) cos(theta_2) 0;...
        0 0 0 1];
    A23 = [...
        1 0 0 0;...
        0 1 0 0;...
        0 0 1 d_3;...
        0 0 0 1];
    
end
A02 = A01*A12;
A0t = A02*A23;

% Rotation matrixes
R01 = A01(1:3,1:3);
R12 = A12(1:3,1:3);
R02 = A02(1:3,1:3);
R2t = A23(1:3,1:3);
R0t = A0t(1:3,1:3);

% Origen location world frame
Ori_1 = A01(1:3,4);
Ori_2 = A02(1:3,4);
Ori_t = A0t(1:3,4);

% World frame vectors
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

% Center of mass in linck frame
r1 = H/2;
r2 = L/2;
r3 = l/2;
rc1 = r1*e3;
rc2 = r2*e1;
rc3 = (q3-r3)*e3;
re1 = H*e3;
re2 = L*e1;
re3 = q3*e3;

% Center of mass in world frame
Rc1 = rc1;
Rc2 = R01*rc2+Ori_1;
Rc3 = R02*rc3+Ori_2;
Rct = Ori_t;

% linear velocite
Rc1_dot = jacobian(Rc1,q)*q_dot.';
Rc2_dot =  jacobian(Rc2,q)*q_dot.';
Rc3_dot =  jacobian(Rc3,q)*q_dot.';
Rct_dot =  jacobian(Rct,q)*q_dot.';
% moment of inertia in link frame
I1 = (m1*H^2)/12;
I2 = (m2*L^2)/12;
I3 = (m3*l^2)/12;
% I3 = (m3*l^2)/12+M*q3^2;
% I_t = M*q3^2;
I_tilde_1 = I1*[1,0,0;
    0,1,0;
    0,0,0];
I_tilde_2 = I2*[0,0,0;
    0,1,0;
    0,0,1];
I_tilde_3 = I3*[1,0,0;
    0,1,0;
    0,0,0];

% Forward Redursion
omega0 = e1*0;
a0 = e1*0;
ac0 = e1*0;
ae0 = e1*0;
% first link
u1 = e3;
omega1 = R01'*omega0+u1*q_dot_1;
a1 = R01'*a0+u1*q_ddot_1+cross(omega1,u1*q_dot_1);
ac1 = R01'*ae0+cross(a1,rc1)+cross(omega1,cross(omega1,rc1));
ae1 = R01'*ae0+cross(a1,re1)+cross(omega1,cross(omega1,re1));
% Second link
u2 = e1;
omega2 = R12'*omega1+u2*q_dot_2;
a2 = R12'*a1+u2*q_ddot_2+cross(omega2,u2*q_dot_2);
ac2 = R12'*ae1+cross(a2,rc2)+cross(omega2,cross(omega2,rc2));
ae2 = R12'*ae1+cross(a2,re2)+cross(omega2,cross(omega2,re2));
% Thered link
u3 = e3;
omega3 = R2t'*omega2+u3*q_dot_3;
a3 = R2t'*a2;
ac3 = R2t'*ae2+cross(a3,re3)+cross(omega3,cross(omega3,rc3))+u3*q_ddot_3+2*cross(omega3,u3)*q_dot_3;
ae3 = R2t'*ae2 + cross(a3,re3)+cross(omega3,cross(omega3,re3))+u3*q_ddot_3+2*cross(omega3,u3)*q_dot_3;

% Backward recursion
% Force on link 3
Fext = -M*g*e3;
Mext = 0*e1;
f3 = m3*ac3-Fext-m3*R0t'*(-g*e3);
M3 = -Mext+cross(rc3,f3)+cross(re3-rc3,-Fext)+I_tilde_3*a3+cross(omega3,I_tilde_3*omega3);
% Force on link 2
f2 = m2*ac2+R2t'*f3-m2*(R02'*e3*(-g));
M2 = R2t*M3+cross(rc2,f2)+cross(re2-rc2,R2t*f3)+I_tilde_2*a2+cross(omega2,I_tilde_2*omega2);
% Force on link 1
f1 = m1*ac1+R12'*f2-m1*R01'*(-g)*e3;
M1 = R12*M2+cross(rc1,f1)+cross(re1-rc1,R12*f2)+I_tilde_1*a1+cross(omega1,I_tilde_1*omega1);
%% Q3 numaric
clear all
global H L l elbow1 elbow2 elbow3
H = 0.2;
L = 0.1;
l = 0.6;
elbow1 = 1;
elbow2 = 1;
elbow3 = 1;
% Calulating mass
d = 0.015;
rou = 7800;
m1 = pi*(d/2)^2*H*rou;
m2 = pi*(d/2)^2*L*rou;
m3 = pi*(d/2)^2*l*rou;
M = 0.5;
g = 9.81;
t = 0:0.001:2;
prof = 3;
% World frame vectors
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];
Force = [];
Moment = [];
for index = t
    q=q_plan(prof,index);
    q_dot=q_dot_plan(prof,index);
    q_dot2=q_dot2_plan(prof,index);
    theta_1 = q(1);
    theta_2 = q(2);
    d_3 = q(3);
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q_dot_1 = q_dot(1);
    q_dot_2 = q_dot(2);
    q_dot_3 = q_dot(3);
    q_ddot_1 = q_dot2(1);
    q_ddot_2 = q_dot2(2);
    q_ddot_3 = q_dot2(3);
    A01 = [...
        cosd(theta_1) -sind(theta_1) 0 0;...
        sind(theta_1) cosd(theta_1) 0 0;...
        0 0 1 H;...
        0 0 0 1];
    A12 = [...
        1 0 0 L;...
        0 cosd(theta_2) -sind(theta_2) 0;...
        0 sind(theta_2) cosd(theta_2) 0;...
        0 0 0 1];
    A23 = [...
        1 0 0 0;...
        0 1 0 0;...
        0 0 1 d_3;...
        0 0 0 1];
    A02 = A01*A12;
    A0t = A02*A23;
    % Rotation matrixes
    R01 = A01(1:3,1:3);
    R12 = A12(1:3,1:3);
    R02 = A02(1:3,1:3);
    R2t = A23(1:3,1:3);
    R0t = A0t(1:3,1:3);
    r1 = H/2;
    r2 = L/2;
    r3 = l/2;
    rc1 = r1*e3;
    rc2 = r2*e1;
    rc3 = (q3-r3)*e3;
    re1 = H*e3;
    re2 = L*e1;
    re3 = q3*e3;
    I1 = (m1*H^2)/12;
    I2 = (m2*L^2)/12;
    I3 = (m3*l^2)/12;
    % I3 = (m3*l^2)/12+M*q3^2;
    % I_t = M*q3^2;
    I_tilde_1 = I1*[1,0,0;
        0,1,0;
        0,0,0];
    I_tilde_2 = I2*[0,0,0;
        0,1,0;
        0,0,1];
    I_tilde_3 = I3*[1,0,0;
        0,1,0;
        0,0,0];
    % Forward Redursion
    omega0 = e1*0;
    a0 = e1*0;
    ac0 = e1*0;
    ae0 = e1*0;
    % first link
    u1 = e3;
    omega1 = R01'*omega0+u1*q_dot_1;
    a1 = R01'*a0+u1*q_ddot_1+cross(omega1,u1*q_dot_1);
    ac1 = R01'*ae0+cross(a1,rc1)+cross(omega1,cross(omega1,rc1));
    ae1 = R01'*ae0+cross(a1,re1)+cross(omega1,cross(omega1,re1));
    % Second link
    u2 = e1;
    omega2 = R12'*omega1+u2*q_dot_2;
    a2 = R12'*a1+u2*q_ddot_2+cross(omega2,u2*q_dot_2);
    ac2 = R12'*ae1+cross(a2,rc2)+cross(omega2,cross(omega2,rc2));
    ae2 = R12'*ae1+cross(a2,re2)+cross(omega2,cross(omega2,re2));
    % Thered link
    u3 = e3;
    omega3 = R2t'*omega2+u3*q_dot_3;
    a3 = R2t'*a2;
    ac3 = R2t'*ae2+cross(a3,re3)+cross(omega3,cross(omega3,rc3))+u3*q_ddot_3+2*cross(omega3,u3)*q_dot_3;
    ae3 = R2t'*ae2 + cross(a3,re3)+cross(omega3,cross(omega3,re3))+u3*q_ddot_3+2*cross(omega3,u3)*q_dot_3;
    
    % Backward recursion
    % Force on link 3
    Fext = -M*g*e3;
    Mext = 0*e1;
    f3 = m3*ac3-Fext-m3*R0t'*(-g*e3);
    M3 = -Mext+cross(rc3,f3)+cross(re3-rc3,-Fext)+I_tilde_3*a3+cross(omega3,I_tilde_3*omega3);
    % Force on link 2
    f2 = m2*ac2+R2t'*f3-m2*(R02'*e3*(-g));
    M2 = R2t*M3+cross(rc2,f2)+cross(re2-rc2,R2t*f3)+I_tilde_2*a2+cross(omega2,I_tilde_2*omega2);
    % Force on link 1
    f1 = m1*ac1+R12'*f2-m1*R01'*(-g)*e3;
    M1 = R12*M2+cross(rc1,f1)+cross(re1-rc1,R12*f2)+I_tilde_1*a1+cross(omega1,I_tilde_1*omega1);
    Force = [Force,f2];
    Moment = [Moment,M2];
end
% Plotting
figure();
h = tiledlayout(3,1);
nexttile
plot(t,Force(1,:))
ylabel('F{_x}[N]')
grid on
set(gca,'FontSize',14);
nexttile
plot(t,Force(2,:))
ylabel('F{_y}[N]')
grid on
set(gca,'FontSize',14);
nexttile
plot(t,Force(3,:))
ylabel('F{_z}[N]')
xlabel('Time[s]')
grid on
set(gca,'FontSize',14);
title(h,'Forces in point b')
set(gca,'FontSize',14);

figure();
h = tiledlayout(3,1);
nexttile
plot(t,Moment(1,:))
ylabel('\tau{_x}[Nm]')
grid on
set(gca,'FontSize',14);
nexttile
plot(t,Moment(2,:))
ylabel('\tau{_y}[Nm]')
grid on
set(gca,'FontSize',14);
nexttile
plot(t,Moment(3,:))
ylabel('\tau{_z}[Nm]')
grid on
xlabel('Time[s]')
grid on
set(gca,'FontSize',14);
title(h,'Moments in point b')
set(gca,'FontSize',14);
%% Q4 Forward dynamics
clear all
tspan = 0:0.0001:2;
global H L l elbow1 elbow2 elbow3 Fe
H = 0.2;
L = 0.1;
l = 0.6;
elbow1 = 1;
elbow2 = 1;
elbow3 = 1;
q_0 = q_plan(3,0);
q = [];
Fe = zeros(6,1);
X_HW4 = [];
for i=tspan
    q = [q ;q_plan(3,i)];
    X_HW4 = [X_HW4; x_plan(3,i)];
end
X0 = [q_0(1);0;q_0(2);0;q_0(3);0];
[t,y] = ode45(@(t, y) state_eq(t, y),tspan,X0);

Joint_traj_plotter("Polynomial", q, [y(:,1) y(:,3) y(:,5)], t)
% Calc Norm and plot
X_num = forward_kin([y(:,1) y(:,3) y(:,5)]);
Err_X2 = (X_num - X_HW4).^2;
X_nrm = sqrt(sum(Err_X2, 2));
full_path_len = 0;
for i=1:length(t)-1
    X_HW4_delta = X_HW4(i+1,:) - X_HW4(i,:);
    full_path_len = full_path_len + sqrt(sum(X_HW4_delta.^2)); % root of sum squared distance between every two points
end
X_nrm_percent = 100 * X_nrm ./ full_path_len;

figure("Name","Q4 Position Error Norm"); 
plot(t, X_nrm_percent); grid on; 
xlabel("$t\ (s)$","Interpreter","latex")
ylabel("$error\ (\%)$","Interpreter","latex")

%% Q5
q_init = inverse_kin([0.4 0 0.81], [1 1 1],0);% changed initial state!!
X0 = [q_init(1);0;q_init(2);0;q_init(3);0]; 
[t,y] = ode45(@(t, y) state_eq(t, y),tspan,X0);

Joint_traj_plotter("Polynomial", q, [y(:,1) y(:,3) y(:,5)], t)

% Calc Norm and plot
X_num = forward_kin([y(:,1) y(:,3) y(:,5)]);
Err_X2 = (X_num - X_HW4).^2;
X_nrm = sqrt(sum(Err_X2, 2));
full_path_len = 0;
for i=1:length(t)-1
    X_HW4_delta = X_HW4(i+1,:) - X_HW4(i,:);
    full_path_len = full_path_len + sqrt(sum(X_HW4_delta.^2)); % root of sum squared distance between every two points
end
X_nrm_percent = 100 * X_nrm ./ full_path_len;

figure("Name","Q5 Position Error Norm"); 
plot(t, X_nrm_percent); grid on; 
xlabel("$t\ (s)$", "Interpreter","latex")
ylabel("$error (\%)$", "Interpreter","latex")
%% Q6
elbow1 = 1;
elbow2 = 1;
elbow3 = 1;
q_init = inverse_kin([0.4 0 0.8], [1 1 1],0); % back to original
X0 = [q_init(1);0;q_init(2);0;q_init(3);0];
Fe(3) = -9.81*0.1;
[t,y] = ode45(@(t, y) state_eq(t, y),tspan,X0); % used new ode func with extra 0.1 effective mass 

Joint_traj_plotter("Polynomial", q, [y(:,1) y(:,3) y(:,5)], t)

% Calc Norm and plot
X_num = forward_kin([y(:,1) y(:,3) y(:,5)]);
Err_X2 = (X_num - X_HW4).^2;
X_nrm = sqrt(sum(Err_X2, 2));
full_path_len = 0;
for i=1:length(t)-1
    X_HW4_delta = X_HW4(i+1,:) - X_HW4(i,:);
    full_path_len = full_path_len + sqrt(sum(X_HW4_delta.^2)); % root of sum squared distance between every two points
end
X_nrm_percent = 100 * X_nrm ./ full_path_len;

figure("Name","Q6 Position Error Norm"); 
plot(t, X_nrm_percent); grid on; 
xlabel("$t\ (s)$", "Interpreter","latex")
ylabel("$error (\%)$", "Interpreter","latex")
