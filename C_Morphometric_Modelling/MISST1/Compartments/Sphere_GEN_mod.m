function E=Sphere_GEN_mod(protocol)
% function [E J] = Sphere_GEN(model, protocol, params_deriv)
% returns the diffusion signal and signal jacobian for the given protocol
% and Sphere tissue model
%
% Substrate: spheres with one radius
% Pulse sequence: Generalized gradient spin echo.
% Signal approximation: Matrix Method
%
% x is the list of model parameters in SI units:
% x(1) is the free diffusivity of the material inside and outside the cylinders.
% x(2) is radius of the sphere
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

if ~isfield(protocol,'GENj')
    protocol.GENj=1;
end

indj=protocol.GENj; % this is the index of the radii. 1 if only one.

tau=protocol.tau;
G=protocol.G;
[M, totsizeG]=size(G);

if protocol.diff==0
    Amat0=protocol.A_sph{indj,1};
    Smat0=protocol.S_sph{indj,1};
    Rmat=protocol.R_sph{indj};
    Amat90=protocol.A_sph{indj,2};
    Smat90=protocol.S_sph{indj,2};
    
elseif protocol.diff==1 % matrices for computing derivatives with respect to diffusivity
    Amat0=protocol.A_sph{indj,1};
    Smat0=protocol.S_sph{indj,1};
    Rmat=protocol.RpertD_sph{indj};
    Amat90=protocol.A_sph{indj,2};
    Smat90=protocol.S_sph{indj,2};
    
elseif protocol.diff==2 % matrices for computing derivatives with respect to radius
    Amat0=protocol.Aperta_sph{indj,1};
    Smat0=protocol.Sperta_sph{indj,1};
    Rmat=protocol.Rperta_sph{indj};
    Amat90=protocol.Aperta_sph{indj,2};
    Smat90=protocol.Sperta_sph{indj,2};
    
else
    error('protocol.diff is not suitable')
end
kDn=protocol.kDn_sph;

E=zeros(M,1);

if isfield(protocol,'smalldel') && isfield(protocol,'delta')
    K = floor((protocol.smalldel+1E-10)./tau)+1;
    dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
    
    %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+G(:,3:3:end).^2);
    Gx = G(:,1:3:end);
    Gy = G(:,2:3:end);
    Gz = G(:,3:3:end);
    cos_theta = zeros(size(Gx));
    sin_phi = zeros(size(cos_theta));
    cos_phi = zeros(size(cos_theta));
    cos_theta(G_mag >0) = Gz(G_mag >0)./G_mag(G_mag >0);
    sin_theta = sqrt(1-cos_theta.^2);
    
    sin_phi(sin_theta ~= 0 & G_mag >0) = Gy(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);
    cos_phi(sin_theta ~= 0 & G_mag >0) = Gx(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);
    
    Gind_mat=round(G_mag./protocol.gstep);
    if protocol.mirror == 1
        for m=1:M
            Smat = real(Smat0) + imag(Smat90)*sin_phi(m,1)+1i*(imag(Smat0)*cos_theta(m,1)+imag(Smat90)*cos_phi(m,1));
            Ecal=Smat';
            for ind=1:K(m)
                Gind=Gind_mat(m,ind);
                if Gind~=0
                    AmatU=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m,ind).*sin_theta(m,ind)+1i*(imag(Amat0)*cos_theta(m,ind)+imag(Amat90)*cos_phi(m,ind)*sin_theta(m,ind));
                    AmatUD=AmatU^Gind;
                    Ecal=Ecal*Rmat*AmatUD;
                elseif Gind==0
                    Ecal=Ecal*Rmat;
                end
            end
            Rmid=Rmat^dstmp(m);
            E(m)=Ecal*Rmid*Ecal';
        end
    elseif protocol.mirror == 0
        for m=1:M
            Smat = real(Smat0) + imag(Smat90)*sin_phi(m,1)+1i*(imag(Smat0)*cos_theta(m,1)+imag(Smat90)*cos_phi(m,1));
            Ecal=Smat';
            ProdMat2 = eye(size(Rmat));
            for ind=1:K(m)
                Gind=Gind_mat(m,ind);
                if Gind~=0
                    AmatU=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m,ind).*sin_theta(m,ind)+1i*(imag(Amat0)*cos_theta(m,ind)+imag(Amat90)*cos_phi(m,ind)*sin_theta(m,ind));
                    AmatUD=AmatU^Gind;
                    AmatUDT = AmatUD';
                    ProdMat2 = ProdMat2*Rmat*AmatUDT; % the matrix product for the second gradient
                    Ecal=Ecal*Rmat*AmatUD;
                elseif Gind==0
                    ProdMat2 = ProdMat2*Rmat;
                    Ecal=Ecal*Rmat;
                end
            end
            Rmid=Rmat^dstmp(m);
            E(m)=Ecal*Rmid*ProdMat2*Smat;
        end
        
    end
else
    K=totsizeG/3;
    % Perpendicular
    
    
    %disp('Angle method using exp(i(k-n)theta)')
    G_mag=sqrt(G(:,1:3:end).^2+G(:,2:3:end).^2+G(:,3:3:end).^2);
    Gx = G(:,1:3:end);
    Gy = G(:,2:3:end);
    Gz = G(:,3:3:end);
    cos_theta = zeros(size(Gx));
    sin_phi = zeros(size(cos_theta));
    cos_phi = zeros(size(cos_theta));
    cos_theta(G_mag >0) = Gz(G_mag >0)./G_mag(G_mag >0);
    sin_theta = sqrt(1-cos_theta.^2);
    
    sin_phi(sin_theta ~= 0 & G_mag >0) = Gy(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);
    cos_phi(sin_theta ~= 0 & G_mag >0) = Gx(sin_theta~=0 & G_mag >0)./G_mag(sin_theta~=0 & G_mag >0)./sin_theta(sin_theta~=0 & G_mag >0);
    
    Gind_mat=round(G_mag./protocol.gstep);
    for m=1:M
        Smat = real(Smat0) + imag(Smat90)*sin_phi(m,1)+1i*(imag(Smat0)*cos_theta(m,1)+imag(Smat90)*cos_phi(m,1));
        Ecal=Smat';
        for ind=1:K
            Gind=Gind_mat(m,ind);
            if Gind~=0
                AmatU=real(Amat0)+(imag(Amat90).*kDn).*sin_phi(m,ind).*sin_theta(m,ind)+1i*(imag(Amat0)*cos_theta(m,ind)+imag(Amat90)*cos_phi(m,ind)*sin_theta(m,ind));
                AmatUD=AmatU^Gind;
                Ecal=Ecal*Rmat*AmatUD;
            elseif Gind==0
                Ecal=Ecal*Rmat;
            end
        end
        E(m)=Ecal*Rmat*Smat;
    end
end
if strcmp(protocol.complex,'complex')
    E=[real(mean(E,2));imag(mean(E,2))];
elseif strcmp(protocol.complex,'real')
    E=real(mean(E,2));
elseif strcmp(protocol.complex,'abs')
    E=abs(mean(E,2));
end
end