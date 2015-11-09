function dydt = sir_ode(t,y,info)
% function dydt = sir_ode(t,y,info)
% 
% SIR model
%
% Original source code for
%
% Modeling post-death transmission of Ebola virus disease (EVD): Challenges for inference and opportunities for control
% Joshua S Weitz and Jonathan Dushoff (in review)
% Preprint available at: arXiv:1411.3435
%
% CC-BY-4.0
dydt=zeros(3,1);
S = y(1);
I = y(2);
R = y(3);
N = sum(y);

dydt(1) = -info.beta*S*I/N;
dydt(2) = info.beta*S*I/N - info.gamma*I;
dydt(3) = info.gamma*I;
