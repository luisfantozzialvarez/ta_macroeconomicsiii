%Declare variables
%Endogenous variables
var x pi i r_n omega;
%Exogenous variables
varexo epsilon_omega epsilon_r;
%Parameters
parameters kappa beta gamma delta_pi delta_x rho_omega rho_r theta nu;

%parameter values
beta = 0.99;
gamma = 1;
nu = 1;
theta = 0.8;
delta_x = 0.25;
delta_pi = 1.5;
rho_omega = 0.5;
rho_r = 0.5;

kappa = ((1-beta*theta)*(1-theta))/(theta*(gamma + nu));

model (linear);
pi = kappa*x + beta*pi(+1); %Phillip's Curve
x = x(+1) -(1/gamma)*(i - pi(+1) - r_n); %IS Curve/Linearised Euler Equation
i = delta_pi*pi + delta_x*x + omega; %Taylor Rule
omega = rho_omega*omega(-1) + epsilon_omega; %Stochastic process for omega
r_n = rho_r*r_n(-1) + epsilon_r; %Stochastic process for natural interest rate
end;

%We set the variance of the shock in order that a 1 sd shock equals a 1
%percentage point shock when computing the IRF, i.e.:
shocks;
var epsilon_omega = 0.01^2;
var epsilon_r = 0.01^2;
end;

stoch_simul(irf = 20);





