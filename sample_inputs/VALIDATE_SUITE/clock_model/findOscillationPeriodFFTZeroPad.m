%INPUTS:
%read in ARdata, stored as: [time, A, R, A(traj2), R(traj2)]
%IMPORTANT: The resolution of the FFT is limited by the max time. so for
%trajs of length 200s, add zero padding, otherwise it will choose between
%22, 25, and 28, so all programs return 25. 
%Zero pad to 5000s gives resolution of ~0.1s at the 25s point. 
     %e.g. zeroPadTime=5000 (in s)
     %OUTPUT:
     %For each A(t) R(t) trajectory, returns the most representative period of oscillations based on the FFT
     %periodMat: matrix of Ntraj x 3
     %each row is for a distinct trajectory. Col1 is period of A. Col2 is period of R. Col3 is lag time between A and R.
     % units of time are based on time vector, e.g. if time is in s, reported periods are in s. 
     %
function [periodMat] = findOscillationPeriodFFTZeroPad(ARdata, zeroPadTime)

[~,n] = size(ARdata);
numsets = (n-1)/2;
periodMat=zeros(numsets,3);
time=ARdata(:,1);
dt=time(2)-time(1);
addtime=[time(end)+dt:dt:zeroPadTime]';%Add time out to ~5000s.
Nzero=length(addtime);
zeropad=zeros(Nzero,1);
totaltime=[time;addtime];
for i=1:numsets
    A=[ARdata(:,2*i);zeropad];
    R=[ARdata(:,2*i+1);zeropad];
	[wlA(i),avA(i), At] = findpeaksfftsub(A,totaltime);
    hold on
    [wlR(i),avR(i), Rt] = findpeaksfftsub(R,totaltime);
    
    lagDiff(i) = lagtime(At(:,2),Rt(:,2),At(:,1));%calculate lag with higher sampled A,R
periodMat(i,:)=[wlA(i), wlR(i), lagDiff(i)];
end

%result = [mean(wlA) mean(wlR) mean(avA) mean(avR) mean(lagDiff)];
%...
 %   std(wlA) std(wlR) std(avA) std(avR) std(lagDiff)];
% disp(ARdata(abs(lagDiff),1));
%
% disp([avA avR]);

% legend('A','R');
% title(['Period A:' num2str(wlA) ', Period R:' num2str(wlR)])
%close all
function [wl,av, At] = findpeaksfftsub(A,distance)

av = mean(A);
r = 10;
y = interp(distance,r); % increase sampling frequency by r
y(y>max(distance)) = [];
y(1) = 0;

A = interp1(distance,A,y);
[ind]=find(isnan(A))
A = A-mean(A);
distance = y;
At=[distance, A];
At(1:10,:)
% Find wavelength of A
L = length(A);
if mod(L,2) ~= 0
    A = A(1:end-1,:);
end

Y = fft(A);
L = length(A);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

difft = diff(distance);
Fs = 1/difft(1);
f = Fs*(0:(L/2))/L;
tvals=1./f;%goes from big to small

arr=find(tvals<100);
display('value where t<100')
stopval=arr(1)
%if you zero pad, the magnitude increases at long times.
[~,i]=max(P1(stopval:end));
index=i+stopval-1
wl = 1/f(index);

%     subplot(1,2,1)
%     plot(distance(1:length(A)),A)
%     xlabel('Length (um)', 'fontsize',14);
%     ylabel('EminDT');
%     subplot(1,2,2)
fnum=20;
figure(fnum);
length(f)

 plot(1./f,P1,'LineWidth',3)
 xlabel('Period (s)', 'fontsize',24);
 ylabel('Magnitude', 'fontsize',24);
     title(['Period ' num2str(wl)])
% pause(0.1)
 set(gcf,'color','w');
 set(gca,'FontSize',24);
 xlim([0 100]);
% axis square

function lagDiff = lagtime(s1,s2,t)

[acor,lag] = xcorr(s2,s1);

[~,I] = max(abs(acor));
lagDiff = lag(I);

dt = diff(t);

Fs = 1/dt(1);
lagDiff = lagDiff/Fs;
