%% GSH for isotropic compression
% given a periodic box in 2D, -L<x,y<L
% apply a small strain rate eps such that the velocity profile is V=(eps/2*x,eps/2*y)
clear
clear variables

% INPUT: strain rate tensor (determined which deformation mode we are studying)
v = .2*[-1 0;0 -1];

%initial values for the time-dependent values
rho0   = .5;
Tg0    = 0;
delta0 = 0;

%solve time-dependent values
[t,y]=ode45(@(t,y)odefun(t,y,v),[0 1],[rho0,Tg0,delta0])
rho   = y(:,1);
Tg    = y(:,2);
delta = y(:,3);
[P_T,P_E]=stress(rho,Tg,delta)

%plot
figure(1)
subplot(2,3,1)
plot(t,rho)
xlabel('t')
ylabel('\rho')
subplot(2,3,2)
plot(t,delta)
xlabel('t')
ylabel('\Delta')
subplot(2,3,3)
plot(t,Tg)
xlabel('t')
ylabel('T_g')
subplot(2,3,4)
plot(t,P_T)
xlabel('t')
ylabel('P_T')
subplot(2,3,5)
plot(t,P_E)
xlabel('t')
ylabel('P_\delta')
subplot(2,3,6)
plot(t,P_T+P_E)
xlabel('t')
ylabel('P')


% solve for time-dependent values
function dy=odefun(t,y,v)
eps = v(1,1)+v(2,2); %volumetric strain rate v_ll
vs2 = sum((v(:)-eps/2*[1;0;0;1]).^2); %deviatoric strain rate squared v_ij*^2
% parameters
RT      =1;
f       =1;
lambda1 =1;
% extract variables for better readability
rho   = y(1);
Tg    = y(2);
delta = y(3); % we need to return the whole elastic strain tensor, not just delta
% ode
dy = [-rho*eps
    -RT*(Tg^2-f*vs2) %dTg
    -eps-lambda1*Tg*delta] %dDelta
end

% get stress assuming no deviatoric strain
function [P_T,P_E]=stress(rho,Tg,delta)
%parameters
rho_b = 0.64;
rho_cp = 0.8;
B_0 = 1;
A_0 = 1;
factor = (max(eps,rho-rho_b)./max(0,rho_cp-rho)); % modified as the factor anthony uses looks wrong
A = A_0*factor;
B = B_0*factor;
g_p = 1; %in reality more complicated
P_T = g_p * Tg.^2;
% get elastic stress
P_E = sqrt(delta).*(B.*delta); % assuming no deviatoric strain
end
