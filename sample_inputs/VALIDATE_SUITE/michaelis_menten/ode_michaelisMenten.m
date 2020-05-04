%Michaelis-Menten
%E+S <->E.S -> E+P
%kf, kb .  kcat
% Concentrations are in units of uM (For E0 and S0)

% time is in s  it is either [T0 Tfinal] or an array of time pts.  

% kf in units of uM(-1)s(-1)
%kb and kcat in s(-1)

function [timepts, conc] = ode_michaelisMenten(E0,S0, kf, kb, kcat, time)


y0=zeros(4,1);
y0(1)=E0;%Enzyme
y0(2)=S0;%substrate
y0(3)=0; %Intermediate: E.S 
y0(4)=0;%Product


%opt=odeset('RelTol',1E-4,'AbsTol',1E-7);
[timepts,conc] = ode23s(@(t,y) odes_mm(t,y,kf, kb,kcat),time,y0);

Eeq=conc(end,1);
Seq=conc(end,2);
Ieq=conc(end,3);
Peq=conc(end,4);


Ptot=Seq+Ieq+Peq;
Etot=Eeq+Ieq;

display('Total substrate before/after') %Ptot should be equal to S0
Ptot
display('Total Enzyme') %Ptot should be equal to E0
Etot



function dy = odes_mm(t,y,kf, kb,kcat)

dy = zeros(4,1);    % a column vector

dy(1) = -kf*y(1)*y(2)+(kb+kcat)*y(3);%E
dy(2) = -kf*y(1)*y(2)+kb*y(3);%S


dy(3) =  +kf*y(1)*y(2)-(kb+kcat)*y(3); %E.S

dy(4) = kcat*y(3); %P




