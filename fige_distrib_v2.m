% Original source code for
%
%Modeling post-death transmission of Ebola virus disease (EVD): Challenges for inference and opportunities for control
%Joshua S Weitz and Jonathan Dushoff (in review)
%Preprint available at: arXiv:1411.3435
%
% CC-BY-4.0
clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'fige_distrib_v2';
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
info.n_E    = 6;        % Number of stages for Gamma
info.T_E    = 11;       % Days exposed, on average
info.b_E    = info.n_E/info.T_E;  % Control parameter for Gamma
x=0:0.01:28;
y=gampdf(x,info.n_E,1/info.b_E);
xdays =1:1:28;
for i=1:length(xdays),  % Cumulative veresion of it
  tmpi=find(x>=(i-1) & x<i);
  ycum(i) = sum(y(tmpi)*0.01);
end
tmph=bar(xdays,ycum);
set(tmph,'facecolor',[0.8 0.8 0.8]);
set(tmph,'edgecolor','k');
set(tmph,'linewidth',2);
set(tmph,'barwidth',0.6);
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
xlim([0 30]);

% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
set(gca,'xtick',[1 5 10 15 20 25]);
set(gca,'xminortick','on');
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

xlabel('Time (days)','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Probability of exposed period','fontsize',20,'verticalalignment','bottom','interpreter','latex');
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
