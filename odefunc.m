function dydt = odefunc(t, y, v, D)
dydt = zeros(7,1);
% Some useful scalars and tensors
v_vol = (1/D)*trace(v)*eye(D);              % volumetric
v_dev = v - v_vol;                       % deviatoric/traceless
v_s = sqrt(sum(sum(v_dev.^2)));                % second invariant
v_ll = trace(v);

% transport coefficients
lambda = 1;
lambda1 = 1;
Rt = 1;
f = 1;
alpha1 = 1;

% evolution variables
rho = y(1);
Tg = y(2);
% u_ij_dev11 = y(3);
% u_ij_dev12 = y(4);
% u_ij_dev21 = y(5);
% u_ij_dev22 = y(6);
u_dev = y(3:6);
u_delta = y(7);                              % delta = -trace(u_ij)


dydt(1) = -rho*v_ll;
dydt(2) = Rt*(-Tg^2 + (f^2)*(v_s^2));
% dydt(3) = v_ij_dev(1,1) - lambda*Tg*u_ij_dev11;
% dydt(4) = v_ij_dev(1,2) - lambda*Tg*u_ij_dev12;
% dydt(5) = v_ij_dev(2,1) - lambda*Tg*u_ij_dev21;
% dydt(6) = v_ij_dev(2,2) - lambda*Tg*u_ij_dev22;
% 
% u_ij_dev  = [u_ij_dev11 u_ij_dev12; u_ij_dev21 u_ij_dev22];

v_dev = reshape(v_dev',4,1);

dydt(3:6) = v_dev - lambda*Tg*u_dev;
dydt(7) = -v_ll + alpha1*sum(sum(u_dev.*v_dev)) - lambda1*Tg*u_delta;
end