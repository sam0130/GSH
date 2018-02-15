%% GSH for isotropic compression
% given a periodic box in 2D, -L<x,y<L
% apply a small strain rate eps such that the velocity profile is V=(eps/2*x,eps/2*y)
clear
close all
eps = -1;
vs = 0;

Tg0 = 0;
delta0 = 0;

[t,y]=ode45(@(t,y) odefun(t,y,vs,eps),[0 1],[Tg0,delta0]);
Tg=y(:,1);
delta=y(:,2);

figure(1)
subplot(1,2,1)
plot(t,delta)
xlabel('t')
ylabel('\Delta')
subplot(1,2,2)
plot(t,Tg)
xlabel('t')
ylabel('T_g')


% solve for Tg, delta
function dy=odefun(t,y,vs,eps)
% parameters
RT=1;
f=1;
lambda1=1;
% extract variables for better readability
Tg=y(1);
delta=y(2);
% ode
dy = [-RT*(Tg^2-f*vs) %dTg
    -eps-lambda1*Tg*delta]; %dDelta
end