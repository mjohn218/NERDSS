%take in a clathrin only simulation, calculate number of bound leg pairs
%return time vs bound legs
function[total]=plot_nfsim_clath_only(traj,linestr,  newPlot)

time=traj(:,1);

freelegs=sum(traj(:,10:12)')';
c0=traj(:,3)*3; %total clathrin *3 is number of legs!
homo=sum(traj(:,4:6)')';
%FOR NFSIM homo dimers are double counted!
hetero=sum(traj(:,7:9)')';

totbound=homo/2+hetero;
n=30

f=figure(n)
if(newPlot=='true')
    ax=axes('Parent',f,'LineWidth',2,'FontSize',30,'XScale','linear')
    hold(ax);
end
plot(time, totbound,linestr,'LineWidth',3);%pairs bound
plot(time, (c0-freelegs)/2,'m--','LineWidth',3);%tot free legs is bound pairs*2

total=[time, totbound];