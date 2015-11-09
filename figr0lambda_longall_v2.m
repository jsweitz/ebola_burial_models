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
tmpfilename = 'figr0lambda_longall_v2';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);
tmppos= [0.2 0.2 0.7 0.7];


set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[460 295 969 511]);

% main data goes here
% Forward simulation of SIR model
info.N=10^6;
info.lambda = 1/28;
info.gammarange = 1./[14 28 42];
info.betarange = info.gammarange.*(1+info.lambda./info.gammarange);
info.R0range = info.betarange./info.gammarange;
info.R0 = 1+info.lambda./info.gammarange;

% Create some exponential "data"
s = RandStream('mt19937ar','Seed',2);
RandStream.setGlobalStream(s);
info.trange = 0:1:6/info.lambda;
% Makes the I(t) be exponential + some noise
% where the noise is approximately normal (but guaranteed to make positive
% case counts)
icaset = exp(info.lambda*info.trange).*exp(randn(1,length(info.trange))*0.2);

% Do the model at early times
tmpcase =0;
tmpc=['r-';'g-';'b-'];
for i=1:length(info.gammarange),
  tmpcase=tmpcase+1;
  tmppos= [0.1 0.1+0.275*(i-1) 0.275 0.275];
  tmpa1 = axes('position',tmppos);
  tmph=plot(info.trange(1:2:end),icaset(1:2:end),'ko');
  set(tmph,'markerfacecolor','k');
  hold on
  info.gamma=info.gammarange(i);
  info.beta = info.betarange(i);
  info.R0 = info.beta/info.gamma;
  options = odeset('reltol',1e-6);
  [t,y]=ode45(@sir_ode,info.trange,[info.N-1 1 0],options,info);
  tmph=plot(t,y(:,2),tmpc(i,:));
  set(tmph,'linewidth',3);
  set(gca,'fontsize',20);
  tmpt = text(25,400,sprintf('${\\cal{R}}_0 = %3.1f$',info.R0));
  set(tmpt,'fontsize',20,'interpreter','latex');
  ylabel('$I(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
  if (tmpcase == 1)
    xlabel('Days, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
    set(gca,'xtick',[0 50 100 150]);
  else
    set(gca,'xtick',[]);
  end
  ylim([0 500]);
  set(gca,'ytick',[0:100:400]);
end

% Now do the long-term model
tmppos= [0.375 0.1 0.5 0.825];
tmpa1 = axes('position',tmppos);
tmpcase =0;
tmpc =['r';'g';'b'];
for i=length(info.gammarange):-1:1,
  tmpcase=tmpcase+1;
  info.gamma=info.gammarange(i);
  info.beta = info.betarange(i);
  info.R0 = info.beta/info.gamma;
  options = odeset('reltol',1e-6);
  [curt,y]=ode45(@sir_ode,[0:1:60/info.lambda],[info.N-1 1 0],options,info);
  % Now show the I(t) signals
  tmpi=find(y(:,2)>0);
  tmph=semilogy(curt(tmpi),y(tmpi,2),tmpc(i,:));
  set(tmph,'linewidth',3);
  hold on
  set(gca,'fontsize',20);
end
xlim([0 1200]);
ylim([10^0 1.5*10^6]);
xlabel('Days, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Infected, $I(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(gca,'yaxislocation','right');
% Plot the data over-top
tmph=semilogy(info.trange(1:5:end),icaset(1:5:end),'ko');
set(tmph,'markerfacecolor','k');


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
tmplh = legend('${\cal{R}}_0 = 2.5$','${\cal{R}}_0 = 2.0$','${\cal{R}}_0 = 1.5$',1);
set(tmplh,'interpreter','latex');
legend('boxoff');

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
