function dydt = odefunc(t, y, v, D)
dydt = zeros(7,1);
% Some useful scalars and tensors
v_vol = (1/D)*trace(v)*eye(D);              % volumetric
v_dev = v - v_vol;                          % deviatoric/traceless
v_s = sqrt(sum(sum(v_dev.^2)));             % second invariant
v_ll = trace(v);

v_dev = reshape(v_dev',4,1);

% transport coefficients
lambda = 1;
lambda1 = 1;
Rt = 1;
f = 1;
alpha1 = 1;

% evolution variables
rho = y(1);
Tg = y(2);
u_dev = y(3:6);
u_delta = y(7);                              % delta = -trace(u_ij)


dydt(1) = -rho*v_ll;
dydt(2) = Rt*(-Tg^2 + (f^2)*(v_s^2));
dydt(3:6) = v_dev - lambda*Tg*u_dev;
dydt(7) = -v_ll + alpha1*sum(sum(u_dev.*v_dev)) - lambda1*Tg*u_delta;
end