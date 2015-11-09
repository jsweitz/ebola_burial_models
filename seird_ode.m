function dydt = seird_ode(t,y,info)
% function dydt = seird_ode(t,y,info)
% 
% SEIR model
% 
% Original source code for
%
% Modeling post-death transmission of Ebola virus disease (EVD): Challenges for inference and opportunities for control
% Joshua S Weitz and Jonathan Dushoff (in review)
% Preprint available at: arXiv:1411.3435
%
% CC-BY-4.0
%
dydt=zeros(5,1);
S = y(1);
E = y(2);
I = y(3);
R = y(4);
D = y(5);
B = y(6);
N = sum(y(1:4));

dydt(1) = -info.beta_I*S*I/N-info.beta_D*S*D/N;
dydt(2) = info.beta_I*S*I/N+info.beta_D*S*D/N - info.sigma*E;
dydt(3) = info.sigma*E - info.gamma*I;
dydt(4) = (1-info.f)*info.gamma*I;
dydt(5) = info.f*info.gamma*I-info.rho*D;
dydt(6)=  info.rho*D;
