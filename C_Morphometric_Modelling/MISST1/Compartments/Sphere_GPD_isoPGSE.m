function [E,J]=Sphere_GPD_isoPGSE(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating a sphere compartment.
% 
% [E,J]=Sphere_GPD_isoPGSE(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of impermeable 
% spheres with a single radius and a diffusion protocol specified in the input
% Substrate: impermeable spheres with a single radius
% Diffusion pulse sequence: isotropic pulsed gradient (isoPGSE) - each 
% gradient waveform has three pairs with orthogonal orientations
% Signal approximation: Gaussian Phase Distribution (GPD)  
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 2 vector of model parameters in SI units for Sphere:
%       x(1) - free diffusivity of the material inside the sphere.
%       x(2) - sphere radius    
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.    
%       protocol.grad_dirs1 - is the gradient direction of the first 
%       gradient pair for each measurement. It has size [N 3] where 
%       N is the number of measurements.
%       protocol.grad_dirs2 - is the gradient direction of the second 
%       gradient pair for each measurement, size [N 3]
%       protocol.grad_dirs3 - is the gradient direction of the third 
%       gradient pair for each measurement, size [N 3]
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - gradient waveform separation, size [1 N]
%       protocol.smalldel - pulse duration of each waveform, size [1 N]
%       protocol.roots_sphere - solutions to the Bessel function equation from 
%       function BesselJ_RootsSphere.m
% x_deriv - a vector of 0s and 1s with the same size as x, indicating
%       which parameters are considered in the Jacobian;
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
% 

roots = protocol.roots_sphere;

% Check the roots array is correct
if(abs(roots(1) -  2.0816)>0.0001)
    error('Looks like the roots array is wrong.  First value should be 2.0816, but is %f', roots(1));
end

dRes=x(1);
R=[x(2)]; 


% get the relevalt pusle sequence parameters from protocol
G = protocol.G';
smalldel = protocol.smalldel'./6; % the duration of each pulse, while smalldel is the total duration of the three gradients
delta = protocol.delta';
grad_dirs1 = protocol.grad_dirs1;
grad_dirs2 = protocol.grad_dirs2;
grad_dirs3 = protocol.grad_dirs3;



% Radial wavenumbers
GAMMA = 2.675987E8; % This is what is used throughout Wuzi.
%GAMMA = 2.6751525E8; % This is the latest best estimate of GAMMA (used in Camino)


l_q=size(grad_dirs1,1);
l_a=numel(R);
k_max=numel(roots);

R_mat=repmat(R,[l_q 1]);
R_mat=R_mat(:);
R_mat=repmat(R_mat,[1 1 k_max]);
R_matSq=R_mat.^2;

root_m=reshape(roots,[1 1 k_max]);
alpha_mat=repmat(root_m,[l_q*l_a 1 1])./R_mat;
amSq=alpha_mat.^2;


%Geometric factor B
Bmat_rep = 2./amSq./(repmat(root_m,[l_q*l_a 1 1]).^2 -2);


delta0mx=repmat(delta,[1,l_a]);
delta0mx_rep = delta0mx(:);
delta0mx_rep = repmat(delta0mx_rep,[1 1 k_max]);

smalldelmx=repmat(smalldel,[1,l_a]);
smalldelmx_rep = smalldelmx(:);
smalldelmx_rep = repmat(smalldelmx_rep,[1 1 k_max]);


Gmx=repmat(G,[1,l_a]);
Gmx_rep = Gmx(:); Gmx_rep = repmat(Gmx_rep ,[1 1 k_max]);


grad_dirs1mx = reshape(grad_dirs1,l_q,1,1,3);
grad_dirs1mx = repmat(grad_dirs1mx,[1,l_a]);
grad_dirs1mx_rep = reshape(grad_dirs1mx,l_q*l_a,1,1,3);
grad_dirs1mx_rep = repmat(grad_dirs1mx_rep ,[1 1 k_max]);


grad_dirs2mx = reshape(grad_dirs2,l_q,1,1,3);
grad_dirs2mx = repmat(grad_dirs2mx,[1,l_a]);
grad_dirs2mx_rep = reshape(grad_dirs2mx,l_q*l_a,1,1,3);
grad_dirs2mx_rep = repmat(grad_dirs2mx_rep ,[1 1 k_max]);


grad_dirs3mx = reshape(grad_dirs3,l_q,1,1,3);
grad_dirs3mx = repmat(grad_dirs3mx,[1,l_a]);
grad_dirs3mx_rep = reshape(grad_dirs3mx,l_q*l_a,1,1,3);
grad_dirs3mx_rep = repmat(grad_dirs3mx_rep ,[1 1 k_max]);

CosPsi12mx_rep = sum(grad_dirs1mx_rep.*grad_dirs2mx_rep,4);
CosPsi13mx_rep = sum(grad_dirs1mx_rep.*grad_dirs3mx_rep,4);
CosPsi23mx_rep = sum(grad_dirs2mx_rep.*grad_dirs3mx_rep,4);

exp_smalldel = exp(-amSq.*dRes.*smalldelmx_rep);
exp_delta0 = exp(-amSq.*dRes.*delta0mx_rep);
exp_delta0_smalldel = exp(-amSq.*dRes.*(delta0mx_rep-smalldelmx_rep));
exp_delta0_2smalldel = exp(-amSq.*dRes.*(delta0mx_rep-2*smalldelmx_rep));
exp_delta0_3smalldel = exp(-amSq.*dRes.*(delta0mx_rep-3*smalldelmx_rep));
exp_delta0_4smalldel = exp(-amSq.*dRes.*(delta0mx_rep-4*smalldelmx_rep));
exp_delta0_5smalldel = exp(-amSq.*dRes.*(delta0mx_rep-5*smalldelmx_rep));
exp_delta0_6smalldel = exp(-amSq.*dRes.*(delta0mx_rep-6*smalldelmx_rep));
am_dRes = amSq.*dRes;
am2_dRes2 = amSq.^2.*dRes.^2;


% Cylindrical restriction
% Sii, gradients with the same orientation
sumGii = 2.*2.*(2.*(am_dRes.*smalldelmx_rep+exp_smalldel-1)./am2_dRes2)...
    +2.*2.*(-(1-exp_smalldel).^2)./am2_dRes2...
    - 4.*((1-exp_smalldel).^2.*exp_delta0_smalldel)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0_2smalldel)./am2_dRes2;
sumterms1 = 3.*Gmx_rep.^2.*sumGii;

% gradient orientation G12
% S13, S14, S23, S19, S110,S29
sumG12 = 2.*4.*((1-exp_smalldel).^2.*exp_smalldel)./am2_dRes2...
    -2.*2.*((1-exp_smalldel).^2.*exp_smalldel.*exp_smalldel)./am2_dRes2...
    -2.*2.*((1-exp_smalldel).^2)./am2_dRes2...
    -4.*((1-exp_smalldel).^2.*exp_delta0.*exp_smalldel)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0.*exp_smalldel.^2)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0)./am2_dRes2;

sumterms1 = sumterms1+ Gmx_rep.*Gmx_rep.*CosPsi12mx_rep.*sumG12;

% gradient orientation 13
%S15, S16, S25, S111, S112, S211

sumG13 = 2.*4.*((1-exp_smalldel).^2.*exp_smalldel.^3)./am2_dRes2...
    - 2.*2.*((1-exp_smalldel).^2.*exp_smalldel.^4)./am2_dRes2...
    - 2.*2.*((1-exp_smalldel).^2.*exp_smalldel.^2)./am2_dRes2...
    - 4.*((1-exp_smalldel).^2.*exp_delta0.*exp_smalldel.^3)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0.*exp_smalldel.^4)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0.*exp_smalldel.^2)./am2_dRes2;

sumterms1 = sumterms1+ Gmx_rep.*Gmx_rep.*CosPsi13mx_rep.*sumG13;


% gradient orientation G23 - sumG12 and G23 are the same
%S35, S36, S45, S311, S312, S411
sumterms1 = sumterms1+ Gmx_rep.*Gmx_rep.*CosPsi23mx_rep.*sumG12;

% gradient orientation G21 (2 from first part, 1 from second part)
sumG21 = -4.*((1-exp_smalldel).^2.*exp_delta0_3smalldel)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0_2smalldel)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0_4smalldel)./am2_dRes2;
sumterms1 = sumterms1+ Gmx_rep.*Gmx_rep.*CosPsi12mx_rep.*sumG21;

% gradient orientation G32 (3 from first part, 2 from second part)
sumterms1 = sumterms1+ Gmx_rep.*Gmx_rep.*CosPsi23mx_rep.*sumG21;

% gradient orientation G31 (3 from first part, 1 from second part)
sumG31 = - 4.*((1-exp_smalldel).^2.*exp_delta0_5smalldel)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0_4smalldel)./am2_dRes2...
    + 2.*((1-exp_smalldel).^2.*exp_delta0_6smalldel)./am2_dRes2;

sumterms1 = sumterms1+ Gmx_rep.*Gmx_rep.*CosPsi13mx_rep.*sumG31;




% % S11
% sumterms1 = 2*(6*Gmx_rep.^2).*(2.*(am_dRes.*smalldelmx_rep+exp_smalldel-1)./am2_dRes2);
% 
% % S12, 34,56 
% 
% for i = 1:6 
%      if  i<=2 
%                 
% %         gd1 = grad_dirs1;
%         gd1 = grad_dirs1mx_rep;
%     elseif i == 3 || i == 4
%       
% %          gd1 = grad_dirs2;
%         gd1 = grad_dirs2mx_rep;
%     elseif i == 5 || i == 6  
%        
% %          gd1 = grad_dirs3;
%         gd1 = grad_dirs3mx_rep;
%     end
%     grad1 = Gmx_rep*(-1).^(i+1);
%     for j = i+1:6     % the same happens for i = 7:12, j = i:12 -> an extra 2 in the formulas  
% %             deltamx=repmat(smalldel*(j-i),[1,l_a]);
% %             deltamx_rep = deltamx(:);
% %             deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);
%             deltamx_rep = smalldelmx_rep.*(j-i);
%             exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
%             
%            
%             
%              if  j<=2              
%                 gd2 = grad_dirs1mx_rep;
%             elseif j == 3 || j == 4 
%               
%                 gd2 = grad_dirs2mx_rep;
%             elseif j == 5 || j == 6 
%                 
%                 gd2 = grad_dirs3mx_rep;
%             end
%             grad2 = Gmx_rep*(-1).^(j+1);
% %             
% %             CosPsi = sum(gd1.*gd2,2);
% %             CosPsimx = repmat(CosPsi,[1,l_a]);
% %             CosPsimx_rep = CosPsimx(:); CosPsimx_rep = repmat(CosPsimx_rep ,[1 1 k_max]);
%             CosPsimx_rep = sum(gd1.*gd2,4);
%             sumterms1 = sumterms1+ 4.*(grad1.*grad2.*CosPsimx_rep).*(exp_delta_smalldel.*(1-exp_smalldel).^2)./am2_dRes2;    
%           
%     end 
%     for j = 7:12
% %          deltamx=repmat(smalldel*(j-i-6)+delta,[1,l_a]);
% %             deltamx_rep = deltamx(:);
% %             deltamx_rep = repmat(deltamx_rep,[1 1 k_max]);
%             deltamx_rep = smalldelmx_rep.*(j-i-6)+delta0mx_rep;
%             exp_delta_smalldel = exp(-amSq.*dRes.*(deltamx_rep-smalldelmx_rep));
%             
%                        
%             if   j ==7 || j == 8
%              
%                 gd2 = grad_dirs1mx_rep;
%             elseif j == 9 || j == 10
%               
%                 gd2 = grad_dirs2mx_rep;
%             elseif  j == 11 || j == 12
%             
%                 gd2 = grad_dirs3mx_rep;
%             end
%             grad2 = -Gmx_rep*(-1).^(j+1); 
%             
% %             CosPsi = sum(gd1.*gd2,2);
% %             CosPsimx = repmat(CosPsi,[1,l_a]);
% %             CosPsimx_rep = CosPsimx(:); CosPsimx_rep = repmat(CosPsimx_rep ,[1 1 k_max]);
%              CosPsimx_rep = sum(gd1.*gd2,4);
%             sumterms1 = sumterms1+ 2.*(grad1.*grad2.*CosPsimx_rep).*(exp_delta_smalldel.*(1-exp_smalldel).^2)./am2_dRes2;  
%          
%     end
% end

 sumterms1 = Bmat_rep.*sumterms1;


testinds = find(sumterms1(:,:,end)>0);
test = sumterms1(testinds,1)./sumterms1(testinds,end);
if(min(test)<1E4)
    warning('Ratio of largest to smallest terms in GPD model sum is <1E4.  May need more terms.');
    x;
end

s = sum(sumterms1,3);
s = reshape(s,[l_q,l_a]);
if(min(s)<0)
    warning('Negative sums found in GPD sum.  Setting to zero.');
    s(find(s<0))=0;
    x;
end
%disp(s.*GmxSq)

logE = -0.5*GAMMA^2.*s;
E_r_matrix = exp(logE);

E_r=sum(E_r_matrix,2);

E=E_r;


% Compute the Jacobian matrix; computed numerically
if(nargout>1)
    dx = 0.00001;
     J = zeros(length(E), length(x));
    if nargin < 3 
        
        for i = 1: length(x); % compute the derivatives for all model parameters
            xpert = x;
            xpert(i) = xpert(i)*(1+dx);
            Epert = Sphere_GPD_isoPGSE(xpert, protocol);
            dEdx = (Epert - E)/(xpert(i)*dx);
            J(:,i) = dEdx;
        end
        
       
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian
    
        for i = 1:length(x_deriv);
            if x_deriv(i) ~= 0                  
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = Sphere_GPD_isoPGSE(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        end
        
       
    end   
end