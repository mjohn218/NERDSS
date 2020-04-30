%Analytical solution for reversible bimolecular association, self and distinct.
%Based on Riccati Equation. 

function[At]=riccati(Atot, Btot, kf, kr, time)
%Atot, Btot units of uM, kf in units of uM-1s-1, kr in units of s-1, time
%in units of s
Keq=kf/kr;
if(Btot==-1)
    display('Self')
    q0=kr*Atot;
    q1=-kr;
    q2=-2*kf;
     Aeq=-1/(4*Keq)+sqrt(1/4/Keq^2+2*Atot/Keq)/2;
else
    display('Distinct')
    q0=kr*Atot;
    q1=-kf*(Btot-Atot)-kr;
    q2=-kf;
    
    b=(Btot-Atot+1/Keq)
    Aeq=-b/2+sqrt(b^2+4*Atot/Keq)/2;

end

R=q1;
   S=q2*q0;
   D=R^2-4*S 
   alpha=R/2+sqrt(R^2-4*S)/2
    beta=R/2-sqrt(R^2-4*S)/2
    
    %alpha=-kr/2+sqrt(kr^2+8*kr*kf*Atot)/2;
    %beta=-kr/2-sqrt(kr^2+8*kr*kf*Atot)/2;
    G=(-beta-Atot*q2)/(Atot*q2+alpha)
    numer=(G*alpha*exp(alpha*time)+beta*exp(beta*time));
    denom=(G*exp(alpha*time)+exp(beta*time));
    At=-numer./denom;
    At=At/q2;
    Aeq
