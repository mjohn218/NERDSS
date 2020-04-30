 % Concentrations are in units of uM (For A0 and M0)
% Length scale s is in units of um
% time is in s  it is either [T0 Tfinal] or an array of time pts.  
% A0(uM)
% k in units of uM(-1)s(-1)
%kb in s(-1)
%VAratio is units of um
%%%%%%%%%%%
%EXAMPLE run:
%[t, outputs]=ode_9species_memlocalize(1, 1, 37.55, 5, 0.001, 0.05, 1, 2, 1, 2, 1, logspace(-3, 1, 200));
%This example returns solutions in conc, in units of uM (same as A0 or M0). 
% OUTPUT:
% timepts (in units of s)
%conc(timepts)=[A, B, M, AM, BM, MAB, ABM, AB, MABM]
% EXAMPLE EQUILIBRIUM:=[0.011579431254713,0.011579431259621, 3.670670706027630,0.085008558198333,0.085008558234364,0.000049217537807, 0.000049217537807, 0.000006704161412, 0.903306871252883]
%%%%%%%%%%%%%%%
function [timepts, conc] = ode_9species_memlocalize(A0,B0, M0, VAratio, s, kfpp, kbpp, kfpm1, kbpm1, kfpm2, kbpm2,time)


y0=zeros(9,1);
y0(1)=A0;%P1
y0(2)=B0;%P2
y0(3)=M0; %M
y0(4)=0;%MP1
y0(5)=0;%P2M
y0(6)=0; %MP1P2
y0(7)=0; %P1P2M
y0(8)=0; %P1P2
y0(9)=0; %MPPM
gamma=VAratio/2/s;% dimensionless!

%opt=odeset('RelTol',1E-4,'AbsTol',1E-7);
[timepts,conc] = ode23s(@(t,y) odes_9species(t,y,kfpp, kbpp, kfpm1, kbpm1,kfpm2, kbpm2, gamma),time,y0);

P1eq=conc(end,1);
P2eq=conc(end,2);
Meq=conc(end,3);
MP1eq=conc(end,4);
P2Meq=conc(end,5);
MP1P2eq=conc(end,6);
P1P2Meq=conc(end,7);
P1P2eq=conc(end,8);
MP1P2Meq=conc(end,9); 

display('Kaeff')
Kaeff=(MP1P2Meq+P1P2Meq+MP1P2eq+P1P2eq)/(P1eq+MP1eq)/(P2eq+P2Meq)
display('Enhancement')
Kaeff/kfpp*kbpp
P1tot=P1eq+P1P2eq+MP1eq+P1P2Meq+MP1P2Meq+MP1P2eq;
P2tot=P2eq+P1P2eq+P2Meq+P1P2Meq+MP1P2Meq+MP1P2eq;

display('Total protein 1') %Ptot should be equal to A0
P1tot
display('Total protein 2') %Ptot should be equal to B0
P2tot
Mtot=Meq+P2Meq+MP1eq+P1P2Meq+MP1P2eq+2*MP1P2Meq;
display('Total Lipids') %Mtot should be equal to M0
Mtot


function dy = odes_9species(t,y,kfpp, kbpp, kfpm1, kbpm1,kfpm2, kbpm2, gamma)

dy = zeros(9,1);    % a column vector

dy(1) = -kfpp*y(1)*y(2)+kbpp*y(8)-kfpp*y(1)*y(5)+kbpp*y(7)-kfpm1*y(1)*y(3)+kbpm1*y(4);%P1
dy(2) = -kfpp*y(1)*y(2)+kbpp*y(8)-kfpp*y(2)*y(4)+kbpp*y(6)-kfpm2*y(2)*y(3)+kbpm2*y(5);%P2


dy(3) =  -kfpm1*y(1)*y(3)+kbpm1*y(4)-kfpm2*y(2)*y(3)+kbpm2*y(5)-gamma*kfpm2*y(3)*y(6)+...
    kbpm2*y(9)-gamma*kfpm1*y(3)*y(7)+kbpm1*y(9)-kfpm1*y(3)*y(8)+kbpm1*y(6)-kfpm2*y(3)*y(8)+kbpm2*y(7); %M

dy(4) = +kfpm1*y(1)*y(3)-kbpm1*y(4)+kbpp*y(9)-gamma*kfpp*y(4)*y(5)-kfpp*y(4)*y(2)+kbpp*y(6); %MP1

dy(5) = +kfpm2*y(2)*y(3)-kbpm2*y(5)+kbpp*y(9)-gamma*kfpp*y(4)*y(5)-kfpp*y(5)*y(1)+kbpp*y(7); %P2M

dy(6) = +kfpm1*y(8)*y(3)-kbpm1*y(6)+kbpm2*y(9)-gamma*kfpm2*y(6)*y(3)+kfpp*y(4)*y(2)-kbpp*y(6); %MP1P2

dy(7) = +kfpm2*y(8)*y(3)-kbpm2*y(7)+kbpm1*y(9)-gamma*kfpm1*y(7)*y(3)+kfpp*y(5)*y(1)-kbpp*y(7); %P1P2M

dy(8) = +kfpp*y(1)*y(2) - kbpp*y(8)+kbpm2*y(7)+kbpm1*y(6)-kfpm1*y(8)*y(3)-kfpm2*y(8)*y(3);%P1P2

dy(9) = -kbpm1*y(9)-kbpm2*y(9)+gamma*kfpm1*y(3)*y(7)+gamma*kfpm2*y(3)*y(6)+gamma*kfpp*y(4)*y(5)-kbpp*y(9); %MP1P2M



