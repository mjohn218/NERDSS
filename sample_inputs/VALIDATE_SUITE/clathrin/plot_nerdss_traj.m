%read in a traj with pip2, ap2, clathrin and synj
%plot different complexes.
%[]=plot_nerdss_traj(traj, dt, newPlot, linestr)
%dt has to be in us
%set newPlot='true' to open new plots, set to 'fals' to add to existing.
%needs to be 4 characters evidently.
%set linestr='r-' or 'b-', to change color


%traj can be an average over multiple runs, just via
%traj=(traj1+traj2+traj3+...trajN)/N;

function[timebound]=plot_nerdss_traj(traj, dt, newPlot, linestr, filename)

%n is the index for the figure+1
n=20
time=traj(:,1)*dt/1E6;
%2,3,4 are free clathrin-cla sites
clfree=sum(traj(:,2:4)')';
%8-13 are bound pairs.
hetero=sum(traj(:,8:10)')';
homo=sum(traj(:,11:13)')';

timebound=[time, hetero+homo];
f=figure(n+1)
maxval=140
if(newPlot=='true')
    ax=axes('Parent',f,'LineWidth',1,'FontSize',12,'XScale','log','TitleFontWeight',...
    'normal','XMinorTick','on','XScale','log', 'YMinorTick','on');
    hold(ax);
    xlim([5E-3 2E2])
    ylim([0 maxval])
     set(ax, 'XTickMode','manual','XTick',[1e-2, 1e0, 1e2],'YTickMode','manual','YTick',[0, maxval/2, maxval],...
        'Units','inches','Position',[0.55 0.55 2.85 1.25]);
end
plot(time, clfree,linestr,'LineWidth',1);
fname=sprintf('unboundCL_%s.eps',filename);
saveas(f,fname,'epsc')
xlabel('time (s)')
ylabel('Bound legs')

f2=figure(n+2)
maxval=140;
if(newPlot=='true')
   ax=axes('Parent',f2,'LineWidth',1,'FontSize',12,'XScale','log','TitleFontWeight',...
    'normal','XMinorTick','on','XScale','log', 'YMinorTick','on');
    hold(ax);
    xlim([5E-3 2E2])
    ylim([0 maxval])
    set(ax, 'XTickMode','manual','XTick',[1e-2, 1e0, 1e2],'YTickMode','manual','YTick',[0, maxval/2, maxval],...
        'Units','inches','Position',[0.55 0.55 2.85 1.25]);
     %'OuterPosition',[0.1 0.1 3.3 1.6],
end

plot(time, hetero,linestr,'LineWidth',1);
plot(time, homo, linestr,'LineWidth',1);
%title('Bound Legs')
fname=sprintf('boundCL_%s.eps',filename);
saveas(f2,fname,'epsc')
xlabel('time (s)')
ylabel('Bound legs')


f3=figure(n+3)
maxval=140;
if(newPlot=='true')
   ax=axes('Parent',f3,'LineWidth',0.5,'FontSize',10,'XScale','log','TitleFontWeight',...
    'normal','XMinorTick','on','XScale','log', 'YMinorTick','on');
    hold(ax);
    xlim([6E-4 2E2])
    ylim([0 maxval])
     set(ax, 'XTickMode','manual','XTick',[1e-2, 1e0, 1e2],'YTickMode','manual','YTick',[0, maxval/2, maxval],...
        'TickDir','out','Units','inches','Position',[0.55 0.55 1.3 0.9],'TickLength',[0.02 0.045]);
end
plot(time, hetero+homo,linestr,'LineWidth',1);
%title('Bound')
fname=sprintf('clacla_%s.eps',filename);
saveas(f3,fname,'epsc')
%xlabel('time (s)')
%ylabel('Bound legs')
