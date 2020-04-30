%READ IN parameter for a 2D bimolecular association A+A-> or A+B->
%Calculate the optimal value for the macroscopic rate constant, 
%kon2D (um^2/s)
%GIVEN:
%ka2D: intrinsic rate (um^2/s)
%sigma: binding radius, molecular length scale in (um)
%Dtot: SUm of diffusion constants of both species (um^2/s)
%SA: Surface area for reaction, a system size is required to define an
%optimal macroscopic rate in 2D!
%NA and NB: copy numbers of each reactant. The length scale is defined
%based on the average separation between NA and NB (or NA for self) in the
%system area.
function [kon2D]=kon2D_value(ka2D, sigma, Dtot, SA, NA, NB)

%length scale b from Area/max(NA,NB)= pi(b^2-s^2), b is the radius of the
%average area per particle
b=2*sqrt(SA/(pi*max(NA,NB))+sigma^2)
bs=b/sigma;
sb=sigma/b;
sb2=sb*sb;

if(b<sigma)
    kon2D=ka2D; %it is unphysical for average separation to be less than particle radius. Default back to maximum rate speed.
else
    ikon=1/ka2D+1/(8*pi*Dtot)* ( -4*log(sb)/(1-sb2)^2-2/(1-sb2)-1);
    kon2D=1.0/ikon;
end
