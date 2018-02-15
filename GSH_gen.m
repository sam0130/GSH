%% GSH for a general compression mode
% given a periodic box in 2D, -L<x,y<L
% General strain rate tensor prescribed: v_ij = [a b; b d]
% Evolution of rho, Tg, u_ij*, \Delta

clear variables; clc; close all;

D = 2;

% INPUT: strain rate tensor
a = -0.2; b = -0.4 ; d = -0.2;
v1 = [a b;
      b d];

v2 = [-0.2 -0.0;    
       -0.0  -0.2 ];   

% v = @(t) v1 + (v2-v1)*heaviside(t-0.5);
   
% Some useful scalars and tensors
% Copy same into odefunc.m
% v_vol = (1/D)*trace(v)*eye(D);                  % volumetric
% v_dev = v - v_vol;                              % deviatoric/traceless
% v_s = sqrt(sum(sum(v_dev.^2)));                 % second invariant
% v_ll = trace(v);

% Initial values of time-dependent variables
rho0 = 0.5;                         % density/volume fraction
Tg0 = 0.5;                            % Granular Temperature
u_dev0 = zeros(D);                  % Deviatoric strain
u_delta0 = 0;                       % Isotropic strain (-u_11 - u_22 - u_33)

y0 = [rho0; Tg0; u_dev0(:); u_delta0]';
tspan1 = [0 1.5];
tspan2 = [1.5 5];

% solve time dependent values
[t1, y1] = ode45(@(t1,y1) odefunc(t1, y1, v1, D), tspan1 , y0);
[t2, y2] = ode45(@(t2,y2) odefunc(t2, y2, v2, D), tspan2 , y1(end,:));
rho = [y1(:,1); y2(:,1)];
Tg = [y1(:,2); y2(:,2)];
u_dev = y(:,3:6);
u_delta = y(:,7);

% Evaluating stresses
% A = 1; B = 1; g_p = 1;  eta_1 = 1;
% u_s = (u_dev(:,1).^2 + u_dev(:,2).^2 + u_dev(:,3).^2 + u_dev(:,4).^2).^0.5; % second invariant of dev strain tensor
% P_delta = (u_delta.^0.5).*(B*u_delta + A*u_s.^2./2./u_delta);               % 1/D*trace(elastic stress)
% stress_elas_vol = P_delta*reshape(eye(2)',1,4);                             % P_delta*I
% stress_elas_dev = (-2*A*u_delta.^0.5).*u_dev;
% stress_elas = stress_elas_vol + stress_elas_dev;
% P_T = (g_p * Tg.^2)*reshape(eye(2)',1,4);                                   % Thermal/fluid pressure
% stress_vis = eta_1*Tg*reshape(v_dev',1,4);                                  % viscous stress
% stress_tot = stress_elas + P_T - stress_vis;                                % total stress tensor

%% plotting
figure
subplot(1,3,1)
plot(t, rho)
xlabel('t')
ylabel('\rho')
subplot(1,3,2)
plot(t, Tg)
xlabel('t')
ylabel('Tg')
subplot(1,3,3)
plot(t, u_delta)
xlabel('t')
ylabel('\Delta')

figure
subplot(2,2,1)
plot(t, u_dev(:,1))
xlabel('t')
ylabel('u_{11}^*')
subplot(2,2,2)
plot(t, u_dev(:,2))
xlabel('t')
ylabel('u_{12}^*')
subplot(2,2,3)
plot(t, u_dev(:,3))
xlabel('t')
ylabel('u_{21}^*')
subplot(2,2,4)
plot(t, u_dev(:,4))
xlabel('t')
ylabel('u_{22}^*')

figure
subplot(2,2,1)
hold on
plot(t, stress_tot(:,1))
plot(t, stress_elas(:,1))
plot(t, P_T(:,1))
plot(t, stress_vis(:,1))
legend('\sigma_{11}','\pi_{11}','P_T','\eta_1T_gv_{11}^*')
hold off
xlabel('t')
%ylabel('\sigma_{11}^*')
subplot(2,2,2)
hold on
plot(t, stress_tot(:,2))
plot(t, stress_elas(:,2))
plot(t, stress_vis(:,2))
legend('\sigma_{12}','\pi_{12}','\eta_1T_gv_{12}^*')
hold off
xlabel('t')
%ylabel('stress')
subplot(2,2,3)
hold on
plot(t, stress_tot(:,3))
plot(t, stress_elas(:,3))
plot(t, stress_vis(:,3))
legend('\sigma_{21}','\pi_{21}','\eta_1T_gv_{21}^*')
hold off
xlabel('t')
%ylabel('\sigma_{21}^*')
subplot(2,2,4)
hold on
plot(t, stress_tot(:,4))
plot(t, stress_elas(:,4))
plot(t, P_T(:,4))
plot(t, stress_vis(:,4))
legend('\sigma_{22}','\pi_{22}','P_T','\eta_1T_gv_{22}^*')
hold off
xlabel('t')
%ylabel('\sigma_{22}^*')

