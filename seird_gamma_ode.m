function dydt = seird_gamma_ode(t,y,info)
% function dydt = seird_gamma_ode(t,y,info)
% 
% SEIR model
% Original source code for
%
% Modeling post-death transmission of Ebola virus disease (EVD): Challenges for inference and opportunities for control
% Joshua S Weitz and Jonathan Dushoff (in review)
% Preprint available at: arXiv:1411.3435
%
% CC-BY-4.0
%
dydt=zeros(5+info.n_E,1);
S = y(1);
En = y(2:info.n_E+1);
I = y(info.n_E+2);
R = y(info.n_E+3);
D = y(info.n_E+4);
B = y(info.n_E+5);
N = sum(y(1:info.n_E+3));  % S E I and R classes

% Model
dydt(1) = -info.beta_I*S*I/N-info.beta_D*S*D/N;
dydt(2) = info.beta_I*S*I/N+info.beta_D*S*D/N - info.b_E*En(1);
for n=2:info.n_E,  % Boxcar method
  dydt(1+n) = info.b_E*En(n-1) - info.b_E*En(n);
end
dydt(info.n_E+2) = info.b_E*En(info.n_E) - info.gamma*I;
dydt(info.n_E+3) = (1-info.f)*info.gamma*I;
dydt(info.n_E+4) = info.f*info.gamma*I-info.rho*D;
dydt(info.n_E+5)=  info.rho*D;
