%read in ARdata, stored as: [time, A, R, A(traj2), R(traj2)]
%Read in a concatenated file where it reports A then R in each column.
%So min size is 2 columns for AR. 
%
%THIS program finds peaks, and then evaluates separations between them,
%uses a threshold for min peak height, and width--so it sometimes makes
%errors!!
%Does not need zero padding or interpolation. 
function[sepsA, sepsR, lags, locA, locR]=calc_peak_sepsAR(ARdata)

%This is set up to ignore peaks that are too narrow (MinPeakWidth)
%It is also set up based on what we know: the peaks are ~20s apart. That
%way it will ignore any small peaks that occur.
%Or can use MinPeakHeight
time=ARdata(:,1);
[len, N]=size(ARdata);
AR=ARdata(:,2:end);
Ntraj=(N-1)/2;
sepsA=[];
sepsR=[];
lags=[];
%   r = 10;
% y = interp(time,r); % increase sampling frequency by r
% y(y>max(time)) = [];
% y(1) = 0;
% newtime=y;
% length(newtime)
for i=1:1:Ntraj
    
    ind=(i-1)*2+1
%interpolation doesn't alter anything.
%A = interp1(time,AR(:,ind),newtime);
%R=interp1(time,AR(:,ind+1),newtime);
A=AR(:,ind);
R=AR(:,ind+1);
newtime=time;
    [pksA,locA]=findpeaks(A,newtime, 'MinPeakWidth',2,'MinPeakHeight',500);
    [pksR,locR]=findpeaks(R,newtime, 'MinPeakWidth',2,'MinPeakHeight',500);
    size(pksA)
    size(pksR)
    %mA(i)=mean(diff(locA));
    %stdA(i)=std(diff(locA));
    sepsA=[sepsA;diff(locA)];
   % medianA(i)=median(diff(locA));
    sepsR=[sepsR;diff(locR)];
    if(length(locR)==length(locA))
        lags=[lags;(locR-locA)];
    else
        display('Num peaks mismatch')
        lags=[lags;(locR-locA(1:end-1))];
    end
end

        