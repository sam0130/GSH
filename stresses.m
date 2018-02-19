function [stress_elas, P_T, stress_vis, stress_tot ] = stresses( v1, y1, D)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Input variables = v_ij, u_delta, u_dev, Tg (simply: v, y1, y2)

global A B g_p eta_1

% State Variables (combined)
rho = y1(:,1);
Tg = y1(:,2);
u_dev = y1(:,3:6) ;
u_delta = y1(:,7);

% Some initial calculations.....
v_vol = (1/D)*trace(v1)*eye(D);              % volumetric
v_dev = v1 - v_vol;                          % deviatoric/traceless

u_s = (u_dev(:,1).^2 + u_dev(:,2).^2 + u_dev(:,3).^2 + u_dev(:,4).^2).^0.5; % second invariant of dev strain tensor

P_delta = (u_delta.^0.5).*(B*u_delta + A*u_s.^2./2./u_delta);               % 1/D*trace(elastic stress)
stress_elas_vol = P_delta*reshape(eye(2)',1,4);                             % P_delta*I
stress_elas_dev = (-2*A*u_delta.^0.5).*u_dev;
stress_elas = stress_elas_vol + stress_elas_dev;
P_T = (g_p * Tg.^2)*reshape(eye(2)',1,4);                                   % Thermal/fluid pressure

stress_vis = eta_1*y1(:,2)*reshape(v_dev',1,4);                             % viscous stress
stress_tot = stress_elas + P_T - stress_vis;                                % total stress tensor

end

