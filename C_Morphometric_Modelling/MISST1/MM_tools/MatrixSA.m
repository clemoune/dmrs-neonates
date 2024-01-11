function [S A kDn]=MatrixSA(q,a,dimA,roo,model)
%
% camino.m--------------------------------------------------------------
% Vector S and matrix A from the Matrix Method formalism (Drobnjak et al, JMR11)
% 
% [S A kDn]=MatrixR(tau,a,dimA,roo,model,D)
% 
% Description: Returns the matrix R which describes the diffusion evolution
% between two gradient point of a diffusion waveform
%
% Parameters: 
% S -  vector which describes the effect of the first gradient point
% A - matrix which describes the effect of a gradient unit (strength gstep, duration tau)
% kDn - a matrix which indicates the sign(m-mprime), as described in (Drobnjak et al, JMR11)
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


% qvec is a projection of the gradient vector onto the x-y plane
if 2*pi*q*a>0.001 
    disp('Warning: Approximation of the exp(i*2*pi*qa) is not accurate enough');
    disp('2*pi*q*a=');
    disp(2*pi*q*a)
    disp('a = ');
    disp(a);
end

kDn=zeros(dimA);
 

if model==0 
    beta = roo(:,5); % normalisation constants
    N = roo(:,1); % index of spherical bessel functions
    M = roo(:,2); % index of legendre polynomials
    alpha = roo(:,4); % roots
elseif model==1
    beta=roo(:,4); 
    N= roo(:,1);
    alpha= roo(:,3); 
elseif model==2
    beta=roo(:,3);
    N=roo(:,1);
    alpha=roo(:,2);
end

if model==0 % SPHERE
    % *******************************
    % Calculate matrix A
    % *******************************
    A0=zeros(dimA,dimA); % for spherical coordinates we need to compute 2 matrices according to the appendix in Drobnjak et al 2011.
    A90=zeros(dimA,dimA);
    S0 = zeros(dimA,1);
    S90 = zeros(dimA,1);
    tol=1e-6;
    % mu are rows of matrices A, S, R and nu are columns.
      if dimA==2 
            alpha2=alpha(2);
            beta1=beta(1);
            beta2=beta(2);
            A0(1,1)=4*pi*beta1^2*a^3/3;
            A90(1,1)=4*pi*beta1^2*a^3/3;
            A0(1,2)=1i*beta1*beta2*q*4*pi^2*(-a^4/alpha2^4)*(3*alpha2*cos(alpha2)+(-3+alpha2^2)*sin(alpha2))*2/3; % calculated in mathematica
            A90(1,2) = 0;
            A0(2,2)=2*pi*beta2^2*a^3/4/alpha2^4*(-2+2*alpha2^2+2*cos(2*alpha2)+alpha2*sin(2*alpha2))*2/3; % calculated in mathematica
            A90(2,2) = A0(2,2);
            A0(2,1)=A0(1,2);
            A90(2,1)=A90(1,2);
            V=4/3*pi*a^3;
            S0(1)=V^(-1/2)*2*pi*beta1*a^3/2*2;
            S0(2)=V^(-1/2)*1i*q*4*pi^2*alpha2*(-a^4/alpha2^4)*(3*alpha2*cos(alpha2)+(-3+alpha2^2)*sin(alpha2))*2/3; % calculated in mathematica
            S90(1)=S0(1);
            S90(2)=0;
            kDn=[0 0;0 0]; % for spherical restricion
      else
          for mu = 1:dimA
              for nu = mu:dimA
                n= N(mu);
                m= M(mu); 
                k=N(nu);
                mprime = M(nu); 
                if m==mprime+1 || m==mprime-1 % follow equation 30 in Drobnjak et al 2011.
                    kDn(mu,nu)=m - mprime;
                    kDn(nu,mu)=mprime-m;
                end
                alpha_mu=alpha(mu);
                alpha_nu=alpha(nu);
                beta_mu=beta(mu);
                beta_nu=beta(nu);
                reA0=0;imA0=0;
                reA90=0;imA90=0;
                if mprime == m 
                    if n == k
                        int_theta = 2/(2*n+1)*factorial(n+m)/factorial(n-m);
                        if n>tol && abs(alpha_mu-alpha_nu)>tol % eigenvalues different and order larger than 0
                            reA0 = 2*pi*beta_mu*beta_nu*(a^3*pi*(alpha_nu*besselj(n-1/2,alpha_nu)*besselj(n+1/2,alpha_mu)-...
                            alpha_mu*besselj(n-1/2,alpha_mu)*besselj(n+1/2,alpha_nu)))/(2*sqrt(alpha_mu*alpha_nu)*(alpha_mu^2-alpha_nu^2))*...
                            int_theta;
                            reA90 = reA0;
                        elseif n>tol && abs(alpha_mu-alpha_nu)<tol && alpha_mu>tol % equal eigenvalues and order larger than 0
                            reA0 = 2*pi*beta_mu^2*a^3*pi/4/alpha_mu*(besselj(n+1/2,alpha_mu)^2-besselj(n-1/2,alpha_mu)*besselj(n+3/2,alpha_mu))*...
                            int_theta;
                            reA90 = reA0;
                        elseif n<tol && abs(alpha_mu-alpha_nu)>tol  % order 0 and different eigenvalues
                            if alpha_mu <tol % only this case necessary
                            reA0 = 2*pi*beta_nu*beta_mu*a^3/alpha_nu^3*(sin(alpha_nu)-alpha_nu*cos(alpha_nu))*...
                            int_theta;
                            else
                            reA0 = 2*pi*beta_nu*beta_mu*a^3/(alpha_mu^3*alpha_nu-alpha_nu^3*alpha_mu)*(alpha_nu*cos(alpha_nu)*sin(alpha_mu)-alpha_mu*cos(alpha_mu)*sin(alpha_nu))*...
                            int_theta;  
                            end
                            reA90 = reA0;
                        elseif n<tol && abs(alpha_mu-alpha_nu)<tol && alpha_mu>tol % order 0 and same eigenvalues
                            reA0 = 2*pi*beta_mu^2*a^3/2/alpha_mu^3*(alpha_mu-cos(alpha_mu)*sin(alpha_mu))*...
                            int_theta;   
                            reA90 = reA0;
                        elseif n<tol && abs(alpha_mu-alpha_nu)<tol && alpha_mu<tol
                            reA0 = 2*pi*beta_mu^2*a^3/3*int_theta;  
                            reA90 = reA0;                            
                        end
                    
                    elseif (n==k+1 ) || (n==k-1 )
                        F1 = @(x) associated_legendre(n,m,x).*associated_legendre(k,m,x).*x;
                        F2 = @(x) spherical_besselj(n,alpha_mu.*x./a).*spherical_besselj(k,alpha_nu.*x./a).*x.^3;
                        imA90 = 0; % due to the cos(theta) term which is 0 or theta = pi/2;
                        imA0 = q*4*pi.^2*beta_mu*beta_nu*integral(F2,0,a)*integral(F1,-1,1);                     
                    end                    
                elseif (m==mprime+1 ) || (m==mprime-1 )
                    if (n==k+1 ) || (n==k-1 )
                        F1 = @(x) associated_legendre(n,m,x).*associated_legendre(k,mprime,x).*sqrt(1-x.^2);
                        F2 = @(x) spherical_besselj(n,alpha_mu.*x./a).*spherical_besselj(k,alpha_nu.*x./a).*x.^3;
                        imA0 = 0;
                        imA90 = q*2*pi^2*beta_mu*beta_nu*integral(F2,0,a)*integral(F1,-1,1);    
                    end
                end   
                A0(mu,nu)=reA0+1i*imA0;
                A90(mu,nu)=reA90+1i*imA90;
                A0(nu,mu)=A0(mu,nu);
                A90(nu,mu)=A90(mu,nu);
              end
          end
      end
      % Calculate S
        for mu = 1:dimA
            reS0 = 0; reS90 = 0; imS0 = 0; imS90 = 0;
            n= N(mu);
            m = M(mu);
            alpha_mu=alpha(mu);
            beta_mu=beta(mu);
            V=4/3*pi*a^3;
            if m == 0
                if n == 0 && alpha_mu <tol
                    reS0 = V^(-1/2)*2*pi*beta_mu*a^3/3*2;
                    reS90 = reS0;
                elseif n == 0 && alpha_mu >tol
                    reS0 = V^(-1/2)*2*pi*beta_mu*(-a^3/alpha_mu^3)*(alpha_mu*cos(alpha_mu)-sin(alpha_mu))*2;
                    reS90 = reS0;
                elseif n == 1
                    imS0 = V^(-1/2)*q*4*pi^2*beta_mu*(-a^4/alpha_mu^4)*(3*alpha_mu*cos(alpha_mu)+(-3+alpha_mu^2)*sin(alpha_mu))*2/3;
                    
                end
            elseif m == 1
                if n == 1
                    imS90 =  V^(-1/2)*q*2*pi^2*beta_mu*(-a^4/alpha_mu^4)*(3*alpha_mu*cos(alpha_mu)+(-3+alpha_mu^2)*sin(alpha_mu))*(-4/3);
                end
                
            
            end
            S0(mu) = reS0+1i*imS0;
            S90(mu) = reS90+1i*imS90;

        end  
        A{1} = A0; A{2} = A90;
        S{1} = S0; S{2} = S90;
   
elseif (model==1) % CYLINDER    
      
        % *******************************
        % Calculate matrix A
        % *******************************
        A=zeros(dimA,dimA);
        S=zeros(dimA,1);
        tol=1e-6;
        % mu are rows of matrices A, S, R and nu are columns.
        if dimA==2 
            alpha2=alpha(2);
            beta1=beta(1);
            beta2=beta(2);
            A(1,1)=pi*beta1*beta1*a^2;
            A(1,2)=1i*a^3.*besselj(2,alpha2).*2*pi^2*beta1*beta2*q/alpha2;
            A(2,2)=pi*beta2.*beta2.*a^2./(2*alpha2).*...
                      (alpha2*besselj(1,alpha2).^2-2*besselj(1,alpha2).*...
                      besselj(2,alpha2)+alpha2*besselj(2,alpha2).^2);
            A(2,1)=A(1,2);
            V=pi*a^2;
            S(1)=pi*V^(-1/2)*beta1*a^2;
            F=@(r) besselj(1,alpha2*r/a).*r.^2;
            S(2)=1i*2*pi^2*V^(-1/2)*beta2*q*integral(F,0,a);
            kDn=[0 -1;1 0];
        else
           for mu = 1:dimA
            for nu=mu:dimA
                reA=0; imA=0;
                n= N(mu);
                k=N(nu);
                if n==k+1 || n==k-1 
                    kDn(mu,nu)=k-n;
                    kDn(nu,mu)=n-k;
                end
                alpha_mu=alpha(mu);
                alpha_nu=alpha(nu);
                beta_mu=beta(mu);
                beta_nu=beta(nu);
                if k==n
                  if n>tol && abs(alpha_mu-alpha_nu)>tol 
                    reA=pi*beta_mu.*beta_nu.*a^2./(alpha_mu^2-alpha_nu^2).*...
                      (alpha_nu*besselj(n-1,alpha_nu).*besselj(n,alpha_mu)-...
                      alpha_mu.*besselj(n-1,alpha_mu).*besselj(n,alpha_nu));
                  elseif n>tol && abs(alpha_mu-alpha_nu)<tol && alpha_mu>tol 
                    reA=pi*beta_mu.*beta_nu.*a^2./(2*alpha_mu).*...
                      (alpha_mu*besselj(n,alpha_mu).^2-2*n*besselj(n,alpha_mu).*...
                      besselj(n+1,alpha_mu)+alpha_mu*besselj(n+1,alpha_mu).^2);
                  elseif n>tol && abs(alpha_mu-alpha_nu)<tol && alpha_mu<tol 
                    reA=0;
                  elseif n<tol && abs(alpha_mu-alpha_nu)>tol 
                      reA=2*pi*beta_mu.*beta_nu.*a^2./(alpha_mu^2-alpha_nu^2).*...
                      (alpha_nu*besselj(n-1,alpha_nu).*besselj(n,alpha_mu)-...
                      alpha_mu.*besselj(n-1,alpha_mu).*besselj(n,alpha_nu));
                  elseif n<tol && abs(alpha_mu-alpha_nu)<tol 
                    reA=pi*beta_mu.*beta_nu.*a^2.*...
                      (besselj(n-1,alpha_mu).^2+besselj(n,alpha_mu).^2);
                  end
                else
                  reA=0;
                end
                if (n==k+1 && k~=0 && k~=-1) || (n==k-1 && k~=0 && k~=1)
                  F=@(r) besselj(n,alpha_mu*r/a).*besselj(k,alpha_nu*r/a).*r.^2;
                  imA=pi^2*beta_mu*beta_nu*q*integral(F,0,a);
                elseif (k==0 && n==1) || (k==1 && n==0)
                   F=@(r) besselj(n,alpha_mu*r/a).*besselj(k,alpha_nu*r/a).*r.^2;
                  imA=2*pi^2*beta_mu*beta_nu*q*integral(F,0,a);
                else
                  imA=0;
                end
                A(mu,nu)=reA+1i*imA;
                A(nu,mu)=A(mu,nu);
            end
           end
        end

        for mu = 1:dimA
           
            n= N(mu);
            alpha_mu=alpha(mu);
            beta_mu=beta(mu);
            V=pi*a^2;
            if n==0
                if alpha_mu == 0
                    S(mu) = 1/sqrt(V)*beta_mu*pi*a^2;
                else
                    S(mu) = 1/sqrt(V)*beta_mu*2*pi*a^2/alpha_mu*besselj(1,alpha_mu);
                end
            elseif n == 1
               S(mu) = 1/sqrt(V)*beta_mu*1i*2*pi^2*a^3/alpha_mu*q*besselj(2,alpha_mu);
            else
               S(mu) = 0;
            end

        end  
       
        
    
elseif model==2 %PLANAR
        % *******************************
        % Calculate matrix A
        % *******************************
     % Calculate A using a Taylor expansion
       A=zeros(dimA,dimA);
        for mu = 1:dimA
           
            for nu=mu:dimA
                n= N(mu);
                k=N(nu);
                alpha_mu=alpha(mu);
                alpha_nu=alpha(nu);
                beta_mu=beta(mu);
                beta_nu=beta(nu);
                if alpha_mu==0 && alpha_nu == 0 
                    A(mu,nu) = beta_mu*beta_nu*2*a;
                elseif alpha_mu == 0
                    if k == 0 % cos
                        A(mu,nu) =  beta_mu*beta_nu*2*a*sin(alpha_nu)/alpha_nu ; %- ...
                     %       2*pi^2*q^2*beta_mu*beta_nu*(2*a^3/alpha_nu^3*(2*alpha_nu*cos(alpha_nu)+(alpha_nu^2-2)*sin(alpha_nu)));  
                    elseif k == 1
                        A(mu,nu) = -beta_mu*beta_nu*4*1i*a^2*pi*q*(alpha_nu*cos(alpha_nu)-sin(alpha_nu))/alpha_nu^2;
                    end
                elseif n == 0 && k == 0 % cos*cos
                    if alpha_mu == alpha_nu
                     A(mu,nu) = beta_mu*beta_nu*a*(1+cos(alpha_mu)*sin(alpha_mu)/alpha_mu);% - ... % 0th order
                      %   2*pi^2*q^2*beta_mu*beta_nu*(a^3/12/alpha_mu^3*(4*alpha_mu^3+6*alpha_mu*cos(2*alpha_mu)+(-3+6*alpha_mu^2)*sin(2*alpha_mu)));  % 2nd order
                    else
                    A(mu,nu) = beta_mu*beta_nu*a*(sin(alpha_mu-alpha_nu)/(alpha_mu-alpha_nu)+sin(alpha_mu+alpha_nu)/(alpha_mu+alpha_nu)); % -...
                     %   2*pi^2*q^2*(beta_mu*beta_nu*2*a^3/(alpha_mu^2-alpha_nu^2)^3*(cos(alpha_nu)*(2*(alpha_mu^4-alpha_nu^4)*cos(alpha_mu)+...
                     %   alpha_mu*(alpha_mu^4-6*alpha_nu^2+alpha_nu^4-2*alpha_mu^2*(1+alpha_nu^2))*sin(alpha_mu))+...
                     %   alpha_nu*(-(alpha_mu^4+alpha_nu^2*(alpha_nu^2-2)-2*alpha_mu^2*(alpha_nu^2+3))*cos(alpha_mu)+...
                     %   4*alpha_mu*(alpha_mu^2-alpha_nu^2)*sin(alpha_mu))*sin(alpha_nu)));
                    end
                elseif n == 0 && k == 1 % cos*sin
                    A(mu,nu) = beta_mu*beta_nu*4*1i*a^2*pi*q/(alpha_mu^2-alpha_nu^2)^2*(alpha_mu*sin(alpha_mu)*(-2*alpha_nu*cos(alpha_nu)+...
                        (alpha_mu-alpha_nu)*(alpha_mu+alpha_nu)*sin(alpha_nu))+cos(alpha_mu)*((alpha_mu-alpha_nu)*(alpha_mu+alpha_nu)*alpha_nu*...
                        cos(alpha_nu)+(alpha_mu^2+alpha_nu^2)*sin(alpha_nu)));
                elseif n == 1 && k == 0 % sin*cos
                    A(mu,nu) = beta_nu*beta_mu*4*1i*a^2*pi*q/(alpha_nu^2-alpha_mu^2)^2*(alpha_nu*sin(alpha_nu)*(-2*alpha_mu*cos(alpha_mu)+...
                        (alpha_nu-alpha_mu)*(alpha_nu+alpha_mu)*sin(alpha_mu))+cos(alpha_nu)*((alpha_nu-alpha_mu)*(alpha_nu+alpha_mu)*alpha_mu*...
                        cos(alpha_mu)+(alpha_nu^2+alpha_mu^2)*sin(alpha_mu)));                      
                elseif n == 1 && k == 1 % sin*sin
                      if alpha_mu == alpha_nu
                     A(mu,nu) = beta_mu*beta_nu*a*(1-cos(alpha_mu)*sin(alpha_mu)/alpha_mu); % - ...
                    %     2*pi^2*q^2*beta_mu*beta_nu*(a^3/12/alpha_mu^3*(4*alpha_mu^3-6*alpha_mu*cos(2*alpha_mu)+(3-6*alpha_mu^2)*sin(2*alpha_mu)));  % 2nd order
                      else
                     A(mu,nu) = beta_mu*beta_nu*2*a*(alpha_nu*cos(alpha_nu)*sin(alpha_mu)-alpha_mu*cos(alpha_mu)*sin(alpha_nu))/(alpha_mu^2-alpha_nu^2); % - ...
                    %    2*pi^2*q^2*( beta_mu*beta_nu*2*a^3/(alpha_mu^2-alpha_nu^2)^3* (sin(alpha_mu)*(alpha_nu*...
                    % (alpha_mu^4-2*alpha_nu^2+alpha_nu^4-2*alpha_mu^2*(3+2*alpha_nu^2))*cos(alpha_nu)+2*(alpha_mu^4-alpha_nu^4)*sin(alpha_nu))+...
                    % cos(alpha_mu)*(4*alpha_mu*(alpha_mu^2-alpha_nu^2)*alpha_nu*cos(alpha_nu)-...
                    % alpha_mu*(alpha_mu^4-6*alpha_nu^2+alpha_nu^4-2*alpha_mu^2*(1+alpha_nu^2))*sin(alpha_nu))));
                      end
                end
                A(nu,mu) = A(mu,nu);
            end
              
        end

        % *******************************
        % Calculate matrix S
        % *******************************
        % analytical expressions
        S=zeros(dimA,1);
        for mu=1:dimA
            n= N(mu);
            alpha_mu=alpha(mu);
            beta_mu=beta(mu);
            V=2*a;
            if alpha_mu == 0
                S(mu) = 1/sqrt(V)*beta_mu*sin(2*pi*a*q)/pi/q;                 
            elseif n==0 % cos
               S(mu) = 1/sqrt(V)*beta_mu*2*a/(alpha_mu^2-4*a^2*pi^2*q^2)*...
                   (alpha_mu*cos(2*a*pi*q)*sin(alpha_mu)-2*a*pi*q*cos(alpha_mu)*sin(2*a*pi*q));
            elseif n==1 % sin
               S(mu) = 1/sqrt(V)*beta_mu*2*1i*a/(alpha_mu^2-4*a^2*pi^2*q^2)*...
                   (2*a*pi*q*cos(2*a*pi*q)*sin(alpha_mu)-alpha_mu*cos(alpha_mu)*sin(2*a*pi*q));
            else
                error('unknown function index - must be 0 or 1')
            end
         
        end

else
        error('Model index not specified. 0 for sphere, 1 for cylinder, 2 for planar.')
end


end

