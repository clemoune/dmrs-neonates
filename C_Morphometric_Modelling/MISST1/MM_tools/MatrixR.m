function R=MatrixR(tau,a,dimA,roo,model,D)
%
% camino.m--------------------------------------------------------------
% Matrix R from the Matrix Method formalism (Drobnjak et al, JMR11)
% 
% R=MatrixR(tau,a,dimA,roo,model,D)
% 
% Description: Returns the matrix R which describes the diffusion evolution
% between two gradient point of a diffusion waveform
%
% Parameters: 
% R - Matrix R from the Matrix Method formalism (Drobnjak et al, JMR11)  
% tau - sampling interval of the gradient waveform
% a - size of restriction (radius for sphere and cylinder restriction,
% half the distance between plates for planar restriction)
% dimA - matrix dimension
% roo - roots
% model - 0 = sphere, 1 = cylinder, 2 = parallel planes 
% D - diffusivity
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Ivana Drobnjak (i.drobnjak@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)

%
R=zeros(dimA,dimA);
for m=1:dimA
  if model==2 % parallel planes
     alpha = roo(m,2);
  elseif model == 1 % cylinder
     alpha = roo(m,3);
  elseif model == 0 % spheres
      alpha = roo(m,4);       
  end
  R(m,m)= exp(-alpha^2*tau*D/(a^2));
end