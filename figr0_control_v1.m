% Original source code for
%
%Modeling post-death transmission of Ebola virus disease (EVD): Challenges for inference and opportunities for control
%Joshua S Weitz and Jonathan Dushoff (in review)
%Preprint available at: arXiv:1411.3435
%
% CC-BY-4.0
clf;
clear all
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figr0_control_v1';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');

% main data goes here
info.lambda_range = 1./[14:7:28];
info.lambda = 1/21;  	% Days
info.n_E    = 6; 	% Number of stages for Gamma
info.T_E    = 11;	% Days exposed, on average
info.b_E    = info.n_E/info.T_E;  % Control parameter for Gamma
info.T_I    = 6;	% Number of days infectious
info.b_I    = 1/info.T_I;% Control parameter for I
info.T_D    = 3;	% Number of days until burial
info.b_D    = 1/info.T_D;% Control parameter for D
info.c_D_range = 0:0.01:1;  % Fraction of secondary infections due to D
tmpc = ['r-';'g-';'b-'];
for l =1:length(info.lambda_range),
  z = -info.lambda_range(l);
  M_E = (info.b_E/(info.b_E-z))^info.n_E;
  M_I = info.b_I/(info.b_I-z);
  M_D = info.b_D/(info.b_D-z);

  for i=1:length(info.c_D_range),
    stats.c_D(i) = info.c_D_range(i);
    stats.c_I(i) = 1-stats.c_D(i);
    stats.Mtot(i) = stats.c_I(i)*M_E*M_I+stats.c_D(i)*M_E*M_I*M_D;
    stats.R0(i) = 1/stats.Mtot(i);
  end
  tmph = plot(stats.c_D,stats.R0,tmpc(l,:));
  set(tmph,'linewidth',3);
  hold on
  tmph2 = plot(stats.c_D,stats.c_I.*stats.R0,tmpc(l,1));
  set(tmph2,'linewidth',3,'LineStyle','--');
  tmpt = text(0.65,1.1*stats.R0(75),sprintf('$\\lambda = 1/%2.0f$',1/info.lambda_range(l)));
  set(tmpt,'interpreter','latex','fontsize',18);
end
xlim([0 1]);
ylim([0 5]);


% loglog(,, '');
%
%
% Some helpful plot commands
% tmph=plot(x,y,'ko');
% set(tmph,'markersize',10,'markerfacecolor,'k');
% tmph=plot(x,y,'k-');
% set(tmph,'linewidth',2);

set(gca,'fontsize',20);

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
% tmplh = legend('','','');
% remove box
% set(tmplh,'visible','off')
% legend('boxoff');

xlabel({'Relative fraction, $\rho_D$,';'of post-death transmission'},'fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Estimated ${\cal{R}}_0$ (--) and ${\cal{R}}_e$ (- -)','fontsize',20,'verticalalignment','bottom','interpreter','latex');
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
