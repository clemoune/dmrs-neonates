function E=AstroSticks_GEN_mod(x, protocol)
% function [E J] = AstroSticks_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and AstroSticks tissue model
%
% Substrate: unifromely oriented sticks 
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: -
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
%
% protocol is a structure containing information related to the diffusion
% sequence
% x_deriv is a vector of 0s and 1s indicating which derivatives are
% calculated
%
% author: Andrada Ianus (a.ianus.11@ucl.ac.uk), Ivana Drobnjak (i.drobnjak@ucl.ac.uk), Daniel
% C. Alexander (d.alexander@ucl.ac.uk)
%
% $Id$

% Model parameters
dRes = x(1);

 if ~isfield(protocol,'GENj')
     protocol.GENj=1;
 end

if ~isfield(protocol,'complex')
  protocol.complex='real';
end

tau=protocol.tau;
G=protocol.G;
[M, totsizeG]=size(G);

GAMMA = 2.675987E8;
% Calculating the parallel & perpendicular signal is similar  (parallel plane & cylindrical restriction )

    K=totsizeG/3;
   %disp('Angle method using exp(i(k-n)theta)')
 
    % get Ndir directions on a spehere
    Ndir = 30;
    fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
    ePerp=ones(M,Ndir); % cylinder restriction
    Bval=zeros(M,Ndir);

    for m=1:M
        for n = 1:Ndir

               fib_dir = fibre_dirs(:,n);            

               G_dot_fibre = G(m,1:3:end)*fib_dir(1)+ G(m,2:3:end)*fib_dir(2) + G(m,3:3:end)*fib_dir(3); 

              for ind=1:K
                 Fpar=sum(G_dot_fibre(1:ind)*tau);
                 Bval(m,n)=Bval(m,n)+(Fpar.^2)*tau; 

              end


        end
    end
    Bval=GAMMA^2*Bval;
    ePar=exp(-Bval.*dRes);
    E = ePar.*ePerp;
    %disp(ePerp)
    if strcmp(protocol.complex,'complex')
      E=[real(mean(E,2));imag(mean(E,2))];
    elseif strcmp(protocol.complex,'real')
      E=real(mean(E,2));
    elseif strcmp(protocol.complex,'abs')
      E=abs(mean(E,2));
    end
end