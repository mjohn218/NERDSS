%
% Usage: input V, A, Kpm1, Kpm2, M0, P01, P02, Kpp and s.
%V is in units um^3, A is in units um^2. Keqs (Kpm1, Kpm2, Kpp) are in units of L/mol. 
%sigma: (s) is in units of um. (e.g. 0.001). 
%Concentrations in units of V^-1: mol/L. Includes P01, P02, and lipids M0!! 
%To change lipid [M0] from units of A^-1 (molecules/um^2) to V^-1 (mol/L), use M0=[M0]*A/V*1E15/6.022E23.
%SC (self-consistent) version iterates on final Kaeff for very small improvement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%IMPORTANT FOR SELF BINDING (homodimers: A+A->C):%%%%%%%%
%SET P02 TO A NEGATIVE NUMBER! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output:
%1. Kaeff (using full theory). units of L/mol. enhance=Kaeff/Kpp.
%2. Kaeff_M0 : Solution of Kaeff using Meq=M0, works when lipids in excess.
%3. All species: Concentrations of all components, units mol/L
%  Allspecies=[Meq, P1, P2, P1M, P2M, P1P2sol,P1P2M, MP1P2, MP1P2M]; 
%  FOR SELF: Allspecies=[Meq, P1, P1M, P1P1sol, P1P1M, MP1P1M];
%4. Percent Complexation: 100*Total protein-protein complexes/ Max possible

%Converge version: Iterates to come up with exact solution by varying
%converge, which determines the correct value of Meq. Meq is the source of
%the approximation in the theory and when it is off, the total lipids in
%complexes does not match total initial lipid concentration: Use this
%difference to update to a better prediction for Meq. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Example Inputs: V=1200 um^3, A=767.66um^2, P01=1E-6mol/L, P02=1E-7 mol/L,
%%M0=2.66E-5mol/L (~25000 lipids/um2), Kpm1=1E6M-1, Kpm2=2E5M-1, Kpp=3E6M-1,
%%s=0.001um. Output: Kaeff=1.8875E9, so enhancement=Kaeff/Kpp=629.17.

function[Kaeff, Kaeff_M0, Allspecies, PercentComplexation]=solveKaeff_converge(Vum, Aum, Kpm1,Kpm2, M0, P01, P02, Kpp, s)

if(P02<0) 
    display('Solving for SELF protein partners: HOMODIMER!');
    self_rxn=1; %solving for protein homodimer: A+A->
    Kpm2=Kpm1;%Only one protein-lipid interaction type.
else
    display('Solving for DISTINCT protein partners: HETERODIMER!');
    self_rxn=0; %solving for protein heterodimer A+B->
end

%gamma
display('Gamma: ')
g=Vum/Aum/s/2

%First calculate simplest estimate of Kaeff (Eq. 3), using Meq=M0
enhance_M0=(g*Kpm1*Kpm2*M0^2+(Kpm1+Kpm2)*M0+1)/(1+Kpm1*M0)/(1+Kpm2*M0);
Kaeff_M0=Kpp*enhance_M0;

%Now calculate Meq, (Eq. 4) to provide complete prediction of Kaeff over all
%conditions.

%Calculate unbound lipids in two limits: Kpp=0 gives Meq0 and Kpp=inf gives MeqCoop.
[Meq0, MeqCoop]=calcUnboundLipid(Vum,Aum, Kpm1,Kpm2, M0, P01, P02, s); %Methods, section A2

%Below is the exact solution to Meq0. This uses the MATLAB 'solve' function
%to solve for PM1 and PM2, which result from a coupled quadratic equation.
%Thus it takes a few seconds. 
if(self_rxn==0)
    [Meq0, ignore, ignore2]=solve_Meq0(Kpm1,Kpm2, M0, P01, P02);
end

MeqCoop=real(MeqCoop);

Meq_curr=Meq0;
 
%Calculate enhancement (Eq. 3) using Meq=Meq0.  
enhance0=(g*Kpm1*Kpm2*Meq_curr^2+(Kpm1+Kpm2)*Meq_curr+1)/(1+Kpm1*Meq_curr)/(1+Kpm2*Meq_curr);%The equation for Kaeff/Ka: Eq. 3
    
 %Calculate lambda (Eq. 8) 
 %Requires the total Pro-Pro complexes formed (using Kaeff instead of Kpp)
  %Divided by max possible complexes
 if self_rxn==0
      PP0=1/2*(P01+P02+1./(Kpp*enhance0))-1/2*sqrt( (P01+P02+1./(Kpp*enhance0)).^2-4*P01*P02 );%Protein complexes
      Pmax=min(P01,P02);%Maximum possible protein-protein complexes
 else
       %Change this equation for self. Assumes that Ka=[AA]/[A]^2, and [A0=A+2*[AA]]  
       PP0=1/2*(P01+1./(4*Kpp*enhance0))-1/2*sqrt( (P01+1./(4*Kpp*enhance0)).^2-P01^2 );%Protein complexes 
       Pmax=P01/2;
 end

  lambda=PP0./Pmax

  %Final Meq is predicted from a linear combination of two limits. 
  %Meq0 (weak Kpp) and MeqCoop (strong Kpp). Switch controlled by amount 
  %of protein complexes formed: lambda.
   Meq=(1-lambda)*Meq0+lambda*MeqCoop; %Eq. 4

 
  %Calculate Enhancement that now uses the final prediction of Meq (unbound lipids)
  enhance=(g*Kpm1*Kpm2*Meq.^2+(Kpm1+Kpm2)*Meq+1)./((1+Kpm1*Meq).*(1+Kpm2*Meq)); %Eq. 3, Kaeff/Kp    
  Kaeff=Kpp.*enhance;
  
  
  %IN This Version, we will iterate over values of Lambda that result in
%proper conservation of total Lipids. If Meq is not accurate enough, the
%number of lipids will not sum to M0

   it=0;
maxIt=1000;%This can be large, the calculation within this loop are simple.
didConverge=0;
while(didConverge==0 && it<maxIt)
   

    Meq=(1-lambda)*Meq0+lambda*MeqCoop; %Eq. 4 
    %Calculate Enhancement that now uses the final prediction of Meq (unbound lipids)
    enhance=(g*Kpm1*Kpm2*Meq.^2+(Kpm1+Kpm2)*Meq+1)./((1+Kpm1*Meq).*(1+Kpm2*Meq)); %Eq. 3, Kaeff/Kpp
    Kaeff=Kpp.*enhance;

    
    %Use Kaeff and initial protein concentrations to predict concentrations of
    %all species.

    if self_rxn==0
       [Allspecies, PercentComplexation]=solve_all_species_hetero(Kaeff, g,Meq,  Kpm1,Kpm2, M0, P01, P02, Kpp);

        %Check whether, based on initial predicted Meq, the lipids sum up to
        %M0. This can fail if Meq is not exact. 
        Meq=Allspecies(1); 
        P1=Allspecies(2);
        P2=Allspecies(3);
        P1M=Allspecies(4);
        P2M=Allspecies(5);
        P1P2sol=Allspecies(6);
        P1P2M=Allspecies(7);
        MP1P2=Allspecies(8);
        MP1P2M=Allspecies(9);
        
        tol=1E-10;%Tolerance for total proteins, they shouldbe correct
        %SET tolerance level for total lipids.
        tolM=0.000001;
        
        P1boundtomem=P1M+P1P2M+MP1P2+MP1P2M;
        P2boundtomem=P2M+P1P2M+MP1P2+MP1P2M;

        %All species from their respective complexes.
        p1sum=P1boundtomem+P1+P1P2sol;%All P1 molecules
        p2sum=P2boundtomem+P2+P1P2sol; %All P2 molecules
        msum=MP1P2+P1P2M+Meq+2*MP1P2M+P1M+P2M; %ALl M molecules

        %Because of the order of equations solved, P1 and P2 are generally correct
        %given any Meq. Test anyways.
        if(P01-p1sum>tol)display('P1 does not add up');end
        if(P02-p2sum>tol)display('P2 does not add up');end
        if(P01-p1sum<-tol)display('P1 exceeds P01!');end
        if(P02-p2sum<-tol)display('P2 exceeds P02!');end

        didConverge=0;
        err=(M0-msum)/M0;
        if(err>tolM)
            %display('M does not add up');
            %msum
            %make lambda smaller.
            scaledown=1-err/100;
            lambda=lambda*scaledown;
        elseif(err<-tolM)
            %display('M exceeds M0!');
            %msum
            %make lambda bigger
            scaleup=1+abs(err)/100;
            lambda=lambda*scaleup;
        else
             didConverge=1;
             display('DID CONVERGE');
             it
             err
             msum
             
        end

    else
        %Self reaction
        [Allspecies, PercentComplexation]=solve_all_species_self(Kaeff, g,Meq,  Kpm1, M0, P01,Kpp);
        M1=Allspecies(1);
        P1=Allspecies(2);
        P1M=Allspecies(3);
        P1P1=Allspecies(4); 
        P1P1M=Allspecies(5);
        MP1P1M=Allspecies(6);
        display('P1 Sum')
        p1sum=P1+P1M+2*P1P1+2*P1P1M+2*MP1P1M
        tol=1E-10;
        tolM=0.01;
        if(P01-p1sum>tol)display('P1 does not add up');end
        if(P01-p1sum<-tol)display('P1 exceeds P01!');end
        
        mSum=P1M+P1P1M+MP1P1M*2+M1;
        err=(M0-mSum)/M0;
        if(err>tolM)
            %display('M does not add up');
            
            %make lambda smaller.
            scaledown=1-err/100;
            lambda=lambda*scaledown;
        elseif(err<-tolM)
            %display('M exceeds M0!');
         
            %make lambda bigger
            scaleup=1+abs(err)/100;
            lambda=lambda*scaleup;
        else
             didConverge=1;
             display('DID CONVERGE');
             it
             err
             mSum
             
        end
        
    end %end of if hetero or homo-dimer reaction

    it=it+1;
    %display('New Lambda')
    %lambda
end%End of while loop.
  
    if(it==maxIt)
        display('FAILED TO CONVERGE')
        lambda
        err
        it
        msum
    end
        
    
end %End of function solveKaeff.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subroutine solve all species hetero

function[Allspecies, PercentComplexation]=solve_all_species_hetero(Kaeff, g,Meq,  Kpm1,Kpm2, M0, P01, P02, Kpp)

TotalComplexes=1/2*(P01+P02+1./Kaeff)-1/2*sqrt( (P01+P02+1./Kaeff).^2-4*P01*P02 );%Protein complexes 
    term=g*Kpm1*Kpm2*Meq.^2+Meq*(Kpm1+Kpm2);
    memfrac=term/(1+term);
    solfrac=1-memfrac;

    TotalMemComplexes=TotalComplexes*memfrac; %[MP1P2M]+[MP1P2]+[P1P2M]

    P1P2sol=TotalComplexes*solfrac;%[P1P2] concentration of complexes in solution.

    %Total P1 and P2 that are not in a protein complex, but are both in solution and
    %membrane.
    P1free=P01-TotalComplexes; %Assume proteins in complex at 1:1 ratio (P1:P2)
    P2free=P02-TotalComplexes;

    memfracP1=Kpm1*Meq./(1+Kpm1*Meq);%[P1M]/ (P1M+P1)
    memfracP2=Kpm2*Meq./(1+Kpm2*Meq);%[P2M]/ (P2M+P2)

    P1M=memfracP1*P1free;%[P1M]
    P2M=memfracP2*P2free;%[P2M]
    P1=P1free-P1M; %[P1]
    P2=P2free-P2M; %[P2]

    P1P2M=P1M*P2*Kpp;
    MP1P2=P2M*P1*Kpp;
    MP1P2M=TotalMemComplexes-P1P2M-MP1P2;
    MP1P2Mcheck=P1M*P2M*Kpp*g;%should be same as above calculation.

    %MP1P2M
    %MP1P2Mcheck

    P1boundtomem=P1M+TotalMemComplexes;
    P2boundtomem=P2M+TotalMemComplexes;
    P1totSol=P1P2sol+P1;
    P2totSol=P1P2sol+P2;


    %Store all species
    %Allspecies=[P1boundtomem, P2boundtomem, P1totSol, P2totSol, TotalComplexes, TotalMemComplexes, P1P2sol, Meq, P1, P2, P1M, P2M];
    Allspecies=[Meq, P1, P2, P1M, P2M, P1P2sol,P1P2M, MP1P2, MP1P2M];
    PercentComplexation=100*(P1P2sol+P1P2M+MP1P2+MP1P2M)/min(P01, P02);%Total proteins in complex relative to max possible protein-protein complexes
end %end solve all species hetero
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subroutine solve all species self

function[Allspecies, PercentComplexation]=solve_all_species_self(Kaeff, g,Meq,  Kpm1, M0, P01,Kpp)

%species for self SI Text section 1B.
    %total complexes has different form.
    TotalComplexes=1/2*(P01+1./(4*Kaeff))-1/2*sqrt( (P01+1./(4*Kaeff)).^2-P01^2 );%Protein complexes 
    
    term=g*Kpm1*Kpm1*Meq.^2+Meq*(Kpm1+Kpm1);
    memfrac=term/(1+term);
    solfrac=1-memfrac;

    TotalMemComplexes=TotalComplexes*memfrac; %[MP1P1M]+[MP1P1]

    P1P1sol=TotalComplexes*solfrac;%[P1P1] concentration of complexes in solution.

    %Total P1 that are not in a protein complex, but are both in solution and
    %membrane.
    P1free=P01-2*TotalComplexes; %two proteins in complex  (P1:P1)
   
    memfracP1=Kpm1*Meq./(1+Kpm1*Meq);%[P1M]/ (P1M+P1)
   
    P1M=memfracP1*P1free;%[P1M]
    P1=P1free-P1M; %[P1]
    
    P1P1M=P1M*P1*2*Kpp;
    
    MP1P1M=TotalMemComplexes-P1P1M;
    MP1P1Mcheck=P1M*P1M*Kpp*g;%should be same as above calculation.

    MP1P1M
    MP1P1Mcheck

    P1boundtomem=P1M+TotalMemComplexes;
    P1totSol=P1P1sol+P1;
    Allspecies=[Meq, P1,  P1M,  P1P1sol, P1P1M, MP1P1M];
    PercentComplexation=100*(P1P1sol+P1P1M+MP1P1M)/(P01/2);%Total proteins in complex relative to max possible protein-protein complexes
end %end solve all species self

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subroutine calcUnboundLipid
%V is in units um^3, A is in units um^2. Keqs (Kpm1, Kpm2, Kpp) are in units of L/mol. 
%sigma: (s) is in units of um. (e.g. 0.001). 
%Concentrations in units of V-1: mol/L. Includes P01, P02, and lipids M0!!! 

function[Meq0, MeqCoop]=calcUnboundLipid(Vum,Aum,  Kpm1,Kpm2, M0, P01, P02, s)

if(P02<0) self_rxn=1; %solving for protein homodimer: A+A->
else self_rxn=0; %solving for protein heterodimer A+B->
end

%Calculate free protein fractions to get the best estimate of Kpm_avg.

if self_rxn==0
    Ptot=P01+P02
    P1=1/2*( P01-M0-1/Kpm1+sqrt( (P01-M0-1/Kpm1)^2+4*P01/Kpm1));%If only P1 binds M0 with Kpm1
    P2=1/2*( P02-M0-1/Kpm2+sqrt( (P02-M0-1/Kpm2)^2+4*P02/Kpm2));%If only P2 binds M0 with Kpm2
    P012=min(P01, P02);
else
    %Changed this for self.
    Ptot=P01
    P1=1/2*( P01-M0-1/Kpm1+sqrt( (P01-M0-1/Kpm1)^2+4*P01/Kpm1));%If only P1 binds M0 with Kpm1
    P2=P1;%Just needed to maintain Kpm_avg equation.
    P012=P01/2;
end
%IF ONLY ONE protein binds lipids
if Kpm1==0
    Ptot=P02; %only one protein2 binds lipid
    P1=0;
end
if Kpm2==0
    Ptot=P01; %only protein1 binds lipids
    P2=0;
end

%Estimate for Kpm when Kpm1 not equal to Kpm2. Works ~<1% error.
%When Kpm1=Kpm2, returns to correct limit, Kpm_avg=Kpm1
Kpm_avg=1/(P1+P2)*(P1*Kpm1+P2*Kpm2)

%Meq0 is after free proteins recruit to membranes via Km
Meq0=1/2*(M0-Ptot-1/Kpm_avg + sqrt( (M0-Ptot-1/Kpm_avg)^2+4*M0/Kpm_avg) );

g=Vum/Aum/s/2; %gamma

Kmav=(Kpm1+Kpm2)/2;
%Equation for MeqCoop is a cubic, aL^3+bL^2+cL+d=0, coefficients below:
a=-g*Kpm1*Kpm2;
b=-2*Kmav+g*Kpm1*Kpm2*M0-2*Kpm1*Kpm2*P012*g;
c=2*M0*Kmav-1-2*Kmav*P012;
d=M0;
%p=[a,b,c,d];
%zs=roots(p);%roots of the polynomial

%explicit definition of the real cubic root.
del0=b^2-3*a*c;
del1=2*b^3-9*a*b*c+27*a^2*d;
c1=((del1+sqrt(del1^2-4*del0^3))/2)^(1/3);
z1=-1/3/a*(b+c1+del0/c1); %Meq with proteins P012 that have two lipid sites

%MeqCoop is after max complexes (each complex has two lipid binding sites) recruit to membranes via Km
MeqCoop=z1; %Eq S14

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here, perform additional minor adjustment, small improvement to predicted 
%Kaeff if P01 not = P02.
%If P01 is much different than P02, there will be
%leftover proteins in solution that can then also bind to lipids, and
%reduce further MeqCoop, since we started with P012 proteins. Does not
%apply to self binding.
if self_rxn==0
    
    Pleft=abs(P01-P02);% proteins in solution that haven't targeted lipids yet.
    if(P02>P01)Kl=Kpm2;
    else Kl=Kpm1;
    end

    M0curr=MeqCoop;
    %calculate how many of these extra proteins are free after trying to bind
    %leftover lipids.
    P1=1/2*( Pleft-M0curr -1/Kl+sqrt( (Pleft-M0curr -1/Kl)^2+4*Pleft/Kl));

    Pbound=Pleft-P1;%proteins bound to lipids
    MeqCoopStar=M0curr-Pbound;%Lipids still unbound, if Pbound=0, no change from MeqCoop.

    
   
    MeqCoop=MeqCoopStar;
end

%When g=1 (no membrane induced cooperativity), MeqCoop=Meq0

%IF ONLY one protein partner binds lipids, no cooperativity is possible.
if Kpm1==0
    MeqCoop=Meq0;
end
if Kpm2==0
    MeqCoop=Meq0;
end
end %end of function calcUnboundLipid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subrouting solve_Meq0
%Reads in affinities between proteins (Kpm1 Kpm2), with
%protein conc defined by P01 and P02
%and lipids conc defined by M0
%units of Kpm1 and P01, P02 and M0 must cancel.
%RETURNS: Equilibrium free lipids and numbers of proteins bound to them.

function[Meq0, PM1, PM2, S]=solve_Meq0(Kpm1,Kpm2, M0, P01, P02)

%At Equilibrium:
%PM1/P1=M*KPM1
%PM2/P2=M*KPM2 
%P01=PM1+P1, P02=PM2+P2, M0=M+PM1+PM2

%PM1/KPM1=(P01-PM1)*(M0-PM1-PM2)
%PM2/KPM2=(PO2-PM2)*(M0-PM1-PM2)

%PM1/KPM1-P01*M0+PM1*P01+PM2*P01+PM1*M0-PM1^2-PM1*PM2
%PM2/KPM2-P02*M0+PM1*P02+PM2*P02+PM2*M0-PM2^2-PM1*PM2

%-PM1^2+PM1*(1/KPM1+P01+M0)-PM1*PM2-P01*M0+PM2*P01
%-PM2^2+PM2*(1/KPM2+P02+M0)-PM1*PM2-P02*M0+PM1*P02

%x=PM1, y=PM2;
syms x y
f1=-x^2+x*(1/Kpm1+P01+M0)-x*y-P01*M0+y*P01;
f2=-y^2+y*(1/Kpm2+P02+M0)-x*y-P02*M0+x*P02;
S=solve(f1==0, f2==0);

%Three (complex) solutions to the above pair of equations. Only one is both
%positive, and produces concentrations less than the initial protein conc.
xs=real(double(S.x));
ys=real(double(S.y));
ind=find(xs>0); %Conc must be positive.
tmp=find(xs(ind)<P01); %Conc cannot exceed initial protein conc!
pair=ind(tmp);
PM1=xs(pair);
PM2=ys(pair);


Meq0=M0-PM1-PM2;
end %end of subroutine solve_Meq0
