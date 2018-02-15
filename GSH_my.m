clear variables; clc; close all;

tspan = [0 1];
e_mag =0.2;                      % strain rate magnitude  
v_ij = e_mag*[-3 0; 0  -1];       % strain rate tensor
v_ll = v_ij(1,1) + v_ij(2,2);
delta0 = 0;
Tg0 = 0; 
rho0 = 0.5;

y0 = [rho0  Tg0 delta0];
[t, y] = ode45(@(t,y) odefunc(t, y, v_ij), tspan , y0);

rho = y(:,1);
Tg = y(:,2);
delta = y(:,3);

p_o = 0.04; gamma_p = 0.033;
p_star = p_o * (delta).*(1 - gamma_p*(delta)); % scaled non-dimensional pressure 



%% plotting
figure
plot(t, y(:,1))
xlabel('t')
ylabel('\rho')
figure
plot(t, y(:,2))
xlabel('t')
ylabel('Tg0')
figure
plot(t, y(:,3))
xlabel('t')
ylabel('\Delta')

figure
plot(delta, p_star)
xlabel('\Delta')
ylabel('p^*')