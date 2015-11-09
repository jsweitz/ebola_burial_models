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
tmpfilename = 'figseird_sleone_v4';
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

%
%
% Case fitting
%
%
x=tdfread('ebola_case_data.txt'); 
tmpi=find(~isnan(x.Cases_SierraLeone) & x.Day>=67 & x.Day<=162);
[p,s]=polyfit(x.Day(tmpi)-x.Day(tmpi(end)),log(x.Cases_SierraLeone(tmpi)),1)
lambda_target = p(1);  % Exponential growth rate as fit to the data
estats.p=p;
estats.s=s;
estats.time=x.Day(tmpi(end:-1:1));
estats.cases_g=x.Cases_SierraLeone(tmpi(end:-1:1));
tmpi=find(~isnan(x.Cases_SierraLeone) & x.Day>=67);
estats.time_all=x.Day(tmpi);
estats.cases_all=x.Cases_SierraLeone(tmpi);

%
%
% Finding R0 using generating functions
%
%
info.n_E    = 6;        % Number of stages for Gamma
info.T_E    = 11;       % Days exposed, on average
info.b_E    = info.n_E/info.T_E;  % Control parameter for Gamma
info.T_I    = 6;        % Number of days infectious
info.b_I    = 1/info.T_I;% Control parameter for I
info.T_D    = [2 4 6];        % Number of days until burial
info.b_D    = 1./info.T_D;% Control parameter for D
info.fracR0_dead_range = 0:0.001:1;  % Fraction of secondary infections due to D
info.f=0.7;  % Fraction of mortality
info.sigma=1/info.T_E;
info.gamma=1/info.T_I;
% Now perform the search
info.lambda_target=lambda_target;
z=-lambda_target;
for n=1:length(info.T_D),
  stats(n).T_D=info.T_D(n);
  stats(n).chi=1/info.T_D(n);
  M_E = (info.b_E/(info.b_E-z))^info.n_E;
  M_I = info.b_I/(info.b_I-z);
  M_D = info.b_D(n)/(info.b_D(n)-z);
  for i=1:length(info.fracR0_dead_range),
    info.fracR0_dead(i) = info.fracR0_dead_range(i);
    info.c_I(i) = 1-info.fracR0_dead(i);
    info.Mtot(i) = info.c_I(i)*M_E*M_I+info.fracR0_dead(i)*M_E*M_I*M_D;
    info.R0(i) = 1/info.Mtot(i);
  end
  stats(n).R0=info.R0;
  stats(n).fracR0_dead = info.fracR0_dead;
  stats(n).betaD = info.fracR0_dead.*info.R0/info.f/stats(n).T_D;
  stats(n).betaI = (1-info.fracR0_dead).*info.R0/info.T_I;
  stats(n).R0_seir= 1/(M_E*M_I);
  stats(n).beta_seir = stats(n).R0_seir/info.T_I;
end
tmpcolbeta = {'r';'g';'c'};

% Focus point identification
% for the n = 2 case (T_D = 4)
tmpcasebeta = [0.05 0.15 0.25];
for i=1:length(tmpcasebeta),
  [tmpbetaval(i) tmpidbeta(i)]=min(abs(stats(2).betaI-tmpcasebeta(i)));
end

% Plotting Upper Left
tmppos= [0.1 0.65 0.275 0.275];
tmpa1 = axes('position',tmppos);
tmpc = {'k--';'k-';'k:'};
for n=1:3,
  tmpi=find(stats(n).betaD>0);
  tmph=plot(stats(n).fracR0_dead(tmpi),stats(n).betaI(tmpi),tmpc{n});
  set(tmph,'linewidth',3);
  hold on
end
n=2;
for i=1:3,
    tmph=plot(stats(n).fracR0_dead(tmpidbeta(i)),stats(n).betaI(tmpidbeta(i)),'ko');
    set(tmph,'markerfacecolor',tmpcolbeta{i},'markersize',12);
end
hold on
tmph = plot(0,stats(n).beta_seir,'kd');
set(tmph,'markerfacecolor','k','markersize',12);
ylabel('$\beta_I$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
ylim([-0.02 0.34]);
set(gca,'xtick',[]);
xlim([-0.1 1.1])
tmplh = legend('$T_D = 2$','$T_D=4$','$T_D=6$',1);
set(tmplh,'interpreter','latex');
legend('boxoff');
set(gca,'fontsize',20);

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
  tmpi=find(double(stats(n).betaD)>0);
  tmph=plot(stats(n).fracR0_dead(tmpi),stats(n).betaD(tmpi),tmpc{n});
  set(tmph,'linewidth',3);
  hold on
end
n=2;
for i=1:3,
    tmph=plot(stats(n).fracR0_dead(tmpidbeta(i)),stats(n).betaD(tmpidbeta(i)),'ko');
    set(tmph,'markerfacecolor',tmpcolbeta{i},'markersize',12);
end
hold on
tmph = plot(0,0,'kd');
set(tmph,'markerfacecolor','k','markersize',12);
ylim([-0.1 1.6]);
xlim([-0.1 1.1]);
set(gca,'ytick',[0.0:0.5:2]);
set(gca,'xtick',[]);
ylabel('$\beta_D$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(gca,'fontsize',20);

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
  tmpi=find(stats(n).betaD>0);
  tmph=plot(stats(n).fracR0_dead(tmpi),stats(n).R0(tmpi),tmpc{n});
  set(tmph,'linewidth',3);
  hold on
end
n=2;
for i=1:3,
  tmph=plot(stats(n).fracR0_dead(tmpidbeta(i)),stats(n).R0(tmpidbeta(i)),'ko');
  set(tmph,'markerfacecolor',tmpcolbeta{i},'markersize',12);
end
tmph = plot(0,stats(n).beta_seir/info.gamma,'kd');
xlim([-0.1 1.1]);
ylim([1.4 2.4]);
set(gca,'ytick',[0:0.25:3]);
set(tmph,'markerfacecolor','k','markersize',12);
ylabel('${\cal{R}}_0$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
xlabel('$\rho_D={\cal{R}}_0(\mathrm{dead~})/{\cal{R}}_0$','fontsize',20,'verticalalignment','top','interpreter','latex');

set(gca,'fontsize',20);


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


% Now do the long-term model
tmpcase =0;
tmpc =['r';'g';'c'];
tmptshift = 7;  % Offset the points by 5 to see them all
n=2;
for i=1:length(tmpcasebeta);
  info.beta_I = tmpcasebeta(i);
  info.beta_D = stats(n).betaD(tmpidbeta(i));
  info.chi=stats(n).chi;
  info.N=6.1*10^6;  % Irrelevant at early times, as long as sufficiently high
  options = odeset('reltol',1e-6);
  [curt,y]=ode45(@seird_gamma_ode,[0:1:12/info.lambda_target],[info.N-1 0 0 0 0 0 0 1 0 0 0],options,info);
  % Now show the I(t) signals
  cumcases = sum(y(:,[8:11])');
  tshift=find(cumcases>exp(estats.p(2)));
  % Case counts
  if (i==1) % Background line all looks the same
    tmph=plot(67+curt(tshift)-curt(tshift(1)),cumcases(tshift),'k-');
    set(tmph,'linewidth',3);
    hold on
  end
  tmph=plot(67+curt(tshift+tmptshift*(i-1):3*tmptshift:end)-curt(tshift(1)),cumcases(tshift+tmptshift*(i-1):3*tmptshift:end),'ko');
  set(tmph,'markerfacecolor',tmpc(i,:),'markersize',12);
  hold on
  set(gca,'fontsize',20);
end
% Plot the "data"
tmph=semilogy(estats.time,estats.cases_g,'k^');
set(tmph,'markersize',12);
[tmpc tmpia]=setdiff(estats.time_all,estats.time);
tmph=semilogy(estats.time_all(tmpia),estats.cases_all(tmpia),'k^');
set(tmph,'markersize',12,'markerfacecolor','k');
hold on

% Modify the view
xlim([40 220]);
set(gca,'xtick',[0:50:250]);
xlabel('Days, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Total cases','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(gca,'yaxislocation','right');
ylim([-20 6000]);
tmpt=text(170,700,{'Sierra';'Leone'});
set(tmpt,'interpreter','latex','fontsize',24);


%
%
% Inset of right panel
% Forward simulation of SEIRD model
% against "data" in log-I(t) space
%
%
tmppos= [0.48 0.55 0.25 0.3];
tmpa4 = axes('position',tmppos);

% Now do the long-term model
tmpcase =0;
tmpc =['r';'g';'c'];
tmptshift = 10;  % Offset the points by 5 to see them all
n=2;
for i=1:length(tmpcasebeta);
  info.beta_I = tmpcasebeta(i);
  info.beta_D = stats(n).betaD(tmpidbeta(i));
  info.chi=stats(n).chi;
  info.N=6.1*10^6;  % Irrelevant at early times, as long as sufficiently high
  options = odeset('reltol',1e-6);
  [curt,y]=ode45(@seird_gamma_ode,[0:1:12/info.lambda_target],[info.N-1 0 0 0 0 0 0 1 0 0 0],options,info);
  % Now show the I(t) signals
  cumcases = sum(y(:,[8:11])');
  % Case counts
  tmph=semilogy(67+curt(tshift)-curt(tshift(1)),cumcases(tshift),'k-');
  set(tmph,'linewidth',3);
  set(tmph,'color',tmpc(i));
  hold on
  set(gca,'fontsize',20);
end
% Plot the "data"
tmph=semilogy(estats.time,estats.cases_g,'k^');
set(tmph,'markersize',12);
[tmpc tmpia]=setdiff(estats.time_all,estats.time);
tmph=semilogy(estats.time_all(tmpia),estats.cases_all(tmpia),'k^');
set(tmph,'markersize',12,'markerfacecolor','k');
hold on

% Clean up the axes
xlim([40 220]);
set(gca,'xtick',[0:50:250]);
xlabel('Days, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Total cases','fontsize',20,'verticalalignment','bottom','interpreter','latex');
ylim([10^0.3 10^3.9]);
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
