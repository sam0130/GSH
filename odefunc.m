function dydt = odefunc(t, y, v, D)
dydt = zeros(7,1);
% Some useful scalars and tensors
s_vol = (1/D)*trace(v)*eye(D);              % volumetric
s_dev = v - s_vol;                          % deviatoric/traceless
s_D = sqrt(sum(sum(s_dev.^2)));             % second invariant
v_ll = trace(v);

s_dev = reshape(s_dev',4,1);

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


dydt(1) = -rho*v_ll;                         % Mass bal   
dydt(2) = Rt*(-Tg^2 + (f^2)*(s_D^2));        % Tg   
dydt(3:6) = s_dev - lambda*Tg*u_dev;         % Deviatoric Strain   
dydt(7) = -v_ll + alpha1*sum(sum(u_dev.*s_dev)) - lambda1*Tg*u_delta;   % Volumetric Strain
end