%%%Model of A->P with kf, and P->A with kb
%Solve analytically--linear system. d/dt[A;P]=[-kf, kb; kf, -kb] [A;P].
%eigenvalues lam=0, -(kf+kb), and eigenvectors, [kb;kf] and [1;-1]
%Coefficient for linear combination are c1=(A0+P0)/(kf+kb), and c2=A0-c1kb,
%or c2=c1kf-P0.
function[at,pt]=reversible_unimolecular(A0, P0, kf, kb)

t=logspace(-7, 2,200); %10^-7:10^2, 200 points.
at=(A0+P0)/(kb+kf)*kb*(1-exp(-tvec*(kb+kf)))+A0*exp(-(kb+kf)*tvec);
pt=(A0+P0)/(kb+kf)*kf*(1-exp(-tvec*(kb+kf)))+P0*exp(-(kb+kf)*tvec);
%steady-state value.
Aeq=kb*(A0+P0)/(kf+kb);
