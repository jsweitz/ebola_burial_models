% Original source code for
%
%Modeling post-death transmission of Ebola virus disease (EVD): Challenges for inference and opportunities for control
%Joshua S Weitz and Jonathan Dushoff (in review)
%Preprint available at: arXiv:1411.3435
%
% CC-BY-4.0
clf;
clear all;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figseird_multiple_v5';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
%set(gcf,'Position', [15 413 1266 393]);
set(gcf,'Position',[460 295 969 511]);

% Cases

%
%
% Panel - Upper Left %
%
%
% Identifiability problem with the betas
J=sym('J',[3 3]);
sigma=sym('sigma');
gamma=sym('gamma');
beta_I=sym('beta_I');
beta_D=sym('beta_D');
rho  = sym('rho');
f_death = sym('f_death');
J(1,1)=-sigma;
J(1,2)=beta_I;
J(1,3)=beta_D;
J(2,1)=sigma;
J(2,2)=-gamma;
J(2,3)=0;
J(3,1)=0;
J(3,2)=f_death*gamma;
J(3,3)=-rho;
[v,d]=eig(J);
info.sigma=1/11;
info.gamma=1/6;
info.f=0.70;  % For Liberia
lambda = d(1,1);
lambda = subs(lambda,{sigma,gamma,f_death},[info.sigma info.gamma info.f]);
lambda_target = 1/21;  % log-doubling every 4 weeks
info.lambda = lambda;
info.lambda_target=lambda_target;
betavec = [0.005:0.005:0.5];  % Range of possible rates of burial
tmpcasebeta = [0.05 0.15 0.25];
for i=1:length(tmpcasebeta),
  tmpidbeta(i)=find(abs(betavec-tmpcasebeta(i))<10^-10);
end
tmpcolbeta = {'r';'g';'c'};

% 3 cases, T_D = 2 4 and 6
fp = fopen('figbeta_betafit_v5.mat');
more off
if (fp<0) % Does not exist
  tmprhovals = 1./[2 4 6];
  for n=1:3,
    clear betaD_exact
    tmplambda = subs(lambda,rho,tmprhovals(n));
    stats(n).rho=tmprhovals(n);
    stats(n).betavec=betavec;
    for i=1:length(betavec),
      i
      tmplambda2 = subs(tmplambda,beta_I,betavec(i));
      betaD_exact(i) = solve(tmplambda2==lambda_target,beta_D);
    end 
  stats(n).betaD_exact=betaD_exact;
  tmpR0 = 1/info.gamma*stats(n).betavec+info.f/stats(n).rho*stats(n).betaD_exact;
  stats(n).fracR0_dead = info.f/stats(n).rho*double(stats(n).betaD_exact)./double(tmpR0);
  stats(n).R0=tmpR0;
  save figbeta_betafit_v5 info stats
  end
else  % If already saved
  fclose(fp);
  load figbeta_betafit_v5
end 
more on

% Plotting Upper Left
tmppos= [0.1 0.65 0.275 0.275];
tmpa1 = axes('position',tmppos);
tmpc = {'k--';'k-';'k:'};
for n=1:3,
  tmpi=find(double(stats(n).betaD_exact)>0);
  tmph=plot(stats(n).fracR0_dead(tmpi),stats(n).betavec(tmpi),tmpc{n});
  set(tmph,'linewidth',3);
  hold on
end
n=2;
for i=1:3,
    tmph=plot(stats(n).fracR0_dead(tmpidbeta(i)),stats(n).betavec(tmpidbeta(i)),'ko');
    set(tmph,'markerfacecolor',tmpcolbeta{i},'markersize',12);
end
% Complex form
beta_seir = ((2*info.lambda_target+(info.sigma+info.gamma))^2-(info.sigma+info.gamma)^2)/(4*info.sigma)+info.gamma
% Simpler form
%beta_seir = (1+info.lambda_target/info.sigma)*(1+info.lambda_target/info.gamma)*info.gamma
hold on
tmph = plot(0,beta_seir,'kd');
set(tmph,'markerfacecolor','k','markersize',12);
ylabel('$\beta_I$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
ylim([-0.02 0.35]);
set(gca,'xtick',[]);
xlim([-0.1 1.1]);
tmplh = legend('$T_D = 2$','$T_D=4$','$T_D=6$',1);
set(tmplh,'interpreter','latex');
legend('boxoff');

%
%
% Panel - Middle Left %
%
%
tmppos= [0.1 0.375 0.275 0.275];
tmpa2 = axes('position',tmppos);
% Plotting
tmpc = {'k--';'k-';'k:'};
for n=1:3,
  tmpi=find(double(stats(n).betaD_exact)>0);
  tmph=plot(stats(n).fracR0_dead(tmpi),double(stats(n).betaD_exact(tmpi)),tmpc{n});
  set(tmph,'linewidth',3);
  hold on
end
n=2;
for i=1:3,
    tmph=plot(stats(n).fracR0_dead(tmpidbeta(i)),double(stats(n).betaD_exact(tmpidbeta(i))),'ko');
    set(tmph,'markerfacecolor',tmpcolbeta{i},'markersize',12);
end
beta_seir = ((2*info.lambda_target+(info.sigma+info.gamma))^2-(info.sigma+info.gamma)^2)/(4*info.sigma)+info.gamma
hold on
tmph = plot(0,0,'kd');
set(tmph,'markerfacecolor','k','markersize',12);
ylim([-0.1 1.6]);
xlim([-0.1 1.1]);
set(gca,'ytick',[0:0.5:2]);
set(gca,'xtick',[]);
ylabel('$\beta_D$','fontsize',20,'verticalalignment','bottom','interpreter','latex');

%
%
% Panel - Bottom Left %
%
%
tmppos= [0.1 0.1 0.275 0.275];
tmpa3 = axes('position',tmppos);
% Plotting
tmpc = {'k--';'k-';'k:'};
for n=1:3,
  tmpi=find(double(stats(n).betaD_exact>0));
  tmph=plot(stats(n).fracR0_dead(tmpi),stats(n).R0(tmpi),tmpc{n});
  set(tmph,'linewidth',3);
  hold on
end
n=2;
for i=1:3,
  tmph=plot(stats(n).fracR0_dead(tmpidbeta(i)),stats(n).R0(tmpidbeta(i)),'ko');
  set(tmph,'markerfacecolor',tmpcolbeta{i},'markersize',12);
end
tmph = plot(0,beta_seir/info.gamma,'kd');
xlim([-0.1 1.1]);
ylim([1.7 2.9]);
set(gca,'ytick',[0:0.25:3]);
set(tmph,'markerfacecolor','k','markersize',12);
ylabel('${\cal{R}}_0$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
xlabel('$\rho_D={\cal{R}}_0(\mathrm{dead~})/{\cal{R}}_0$','fontsize',20,'verticalalignment','top','interpreter','latex');


%
%
% Panel - Right
%
%
% main data goes here
% Forward simulation of SEIRD model
% against "data" in I(t) space
tmppos= [0.375 0.1 0.5 0.825];
tmpa4 = axes('position',tmppos);

% Create some exponential "data"
info.trange = 0:1:12/info.lambda_target;
% Makes the I(t) be exponential (to illustrate the rise)
icaset = exp(info.lambda_target*info.trange);
% Plot the "data"
tmph=plot(info.trange(1:5:end),icaset(1:5:end),'k-');
set(tmph,'linewidth',3);
hold on

% Now do the long-term model
tmpcase =0;
tmpc =['r';'g';'c'];
tmptshift = 5;  % Offset the points by 5 to see them all
n=2;
for i=1:length(tmpcasebeta);
  info.beta_I = tmpcasebeta(i);
  info.beta_D = double(stats(n).betaD_exact(tmpidbeta(i)));
  info.rho=stats(n).rho;
  info.N=10^6;
  options = odeset('reltol',1e-6);
  [curt,y]=ode45(@seird_ode,[0:1:60/info.lambda_target],[info.N-1 0 1 0 0 0],options,info);
  % Now show the I(t) signals
  tmpi=find(y(:,3)>0);
  % Case counts
  tshift=find(y(:,3)>1);
  tvals = curt(tmpi)-curt(tshift(1));
  tmph=plot(curt(tmpi(tshift)+tmptshift*(i-1):3*tmptshift:end)-curt(tshift(1)),y(tmpi(tshift)+tmptshift*(i-1):3*tmptshift:end,3),'ko');
  % Total counts
  %  tmph=semilogy(curt(tmpi),(sum(y(tmpi,2:6)'))'/info.N,tmpc(i,:));
  set(tmph,'markerfacecolor',tmpc(i));
  hold on
  set(gca,'fontsize',20);
end
xlim([-10 160]);
set(gca,'xtick',[0:50:250]);
xlabel('Days, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Infected, $I(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(gca,'yaxislocation','right');
ylim([-20 1600]);

%
%
% Inset of right panel
% Forward simulation of SEIRD model
% against "data" in log-I(t) space
%
%
tmppos= [0.48 0.55 0.3 0.3];
tmpa4 = axes('position',tmppos);

% Create some exponential "data"
info.trange = 0:1:12/info.lambda_target;
% Makes the I(t) be exponential (to illustrate the rise)
icaset = exp(info.lambda_target*info.trange);
% Plot the "data"
tmph=semilogy(info.trange(1:5:end),icaset(1:5:end),'k-');
set(tmph,'linewidth',3);
hold on

% Now do the long-term model
tmpcase =0;
tmpc =['r';'g';'c'];
tmptshift = 10;  % Offset the points by 5 to see them all
n=2;
for i=1:length(tmpcasebeta);
  info.beta_I = tmpcasebeta(i);
  info.beta_D = double(stats(n).betaD_exact(tmpidbeta(i)));
  info.rho=stats(n).rho;
  info.N=10^6;
  options = odeset('reltol',1e-6);
  [curt,y]=ode45(@seird_ode,[0:1:12/info.lambda_target],[info.N-1 0 1 0 0 0],options,info);
  % Now show the I(t) signals
  tmpi=find(y(:,3)>0);
  % Case counts
  tshift=find(y(:,3)>1);
  tvals = curt(tmpi)-curt(tshift(1));
  tmph=semilogy(curt(tmpi(tshift)+tmptshift*(i-1):3*tmptshift:end)-curt(tshift(1)),y(tmpi(tshift)+tmptshift*(i-1):3*tmptshift:end,3),'ko');
  % Total counts
  %  tmph=semilogy(curt(tmpi),(sum(y(tmpi,2:6)'))'/info.N,tmpc(i,:));
  set(tmph,'markerfacecolor',tmpc(i));
  hold on
  set(gca,'fontsize',20);
end
xlim([-10 160]);
set(gca,'xtick',[0:50:250]);
xlabel('Days, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Infected, $I(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
ylim([10^-0.5 10^3.8]);
set(gca,'ytick',10.^[0:3]);



% loglog(,, '');
%
%
% Some helpful plot commands
% tmph=plot(x,y,'ko');
% set(tmph,'markersize',10,'markerfacecolor,'k');
% tmph=plot(x,y,'k-');
% set(tmph,'linewidth',2);


% for use with layered plots
% set(gca,'box','off')

% adjust limits
% tmpv = axis;
% axis([]);
% ylim([]);
% xlim([]);

% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
% set(gca,'xtick',[1 100 10^4])
% set(gca,'ytick',[1 100 10^4])

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);

% title('','fontsize',24)
% 'horizontalalignment','left');

% for writing over the top
% coordinates are normalized again to (0,1.0)
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
% first two points are normalized x, y positions
% text(,,'','Fontsize',14);

% automatic creation of postscript
% without name/date
psprintc(tmpfilenoname);
psprint(tmpfilebwname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*
