function [td_synth, td_diff, params]=svdfid(td, td_fraction, bw, lowf, highf, rank,lo_lw, hi_lw)

% td = FID file
% td_fraction = fraction of the FID used for the svd >> total_nb_FID_points/td_fraction (8 is a good number when there are 2048 points)
% bw = bandwidth (5000 Hz for the 11.7T)
% lowf = lower bound for svd (water peak is centered in 0)
% highf = higher bound for svd (water peak is centered in 0)
% rank = maximum number of peaks to detect
% lo_lw = minimum FWHM of peak to detect
% hi_lx = maximum FWHM of peak to detect



% SVD for parametric estimation of spectra.


delta=1/bw;

% bw in radians/sec

npoints=length(td);
npoints_svd=npoints/td_fraction;
low_svd=lowf*2*pi;
high_svd=highf*2*pi;
% creating Toeplitz matrix from the data (step 1)

T=zeros(npoints_svd,npoints_svd);
% U=zeros(npoints_svd,npoints_svd);
% S=zeros(npoints_svd,npoints_svd);
% V=zeros(npoints_svd,npoints_svd);


for ii=1:npoints_svd
    
    T(ii,:)=td(npoints_svd+1-ii:2*npoints_svd-ii,1);
end


% singular value decomposition (step 2)

[U,S,V]=svd(T);

% reducing rank (step 3)

% Tk=zeros(rank,rank);
% Uk=U(:, 1:rank);
% Sk=S(1:rank, 1:rank);

Vk=V(:, 1:rank);

% Tk=Uk*Sk*Vk';

% least square solution (step 4)

VtK=Vk(2:npoints_svd,:);
VbK=Vk(1:npoints_svd-1,:);

EH=lscov(VtK, VbK);

poles=eig(EH');

decayparms=log(poles)/delta;

% finding the relevant decay parameters for the water region

params_ix=find((imag(decayparms))> low_svd & (imag(decayparms) < high_svd) & (real(decayparms)>lo_lw*pi) &(real(decayparms)<hi_lw*pi));
rank_est=length(params_ix);

svd_params=zeros(1,rank_est);

if rank_est > 0
    
    for ii=1:rank_est
        svd_params(ii)=decayparms(params_ix(ii));
    end
    
    
    % solving for the amplitudes and phases (step 5)
    % first, generating the model matrix
    
    Z=zeros(2048, rank_est);
    for ii=1:2048
        Z(ii,:)=exp(-(svd_params)*delta*(ii-1));
    end
    
    % then solving the least squares equation
    
    amps=lsqr(Z(1:npoints,:), td(1:npoints));
    cond_num=cond(Z(1:npoints,:)'*Z(1:npoints,:))
    
    % creating the synthetic time series
    
    td_synth=zeros(2048,1);
    
    for ii=1:rank_est
        
        td_synth=td_synth+amps(ii)*Z(:,ii);
        
    end
    td_diff=td-td_synth(1:npoints);
    
    
else
    td_diff=td;
end


params=zeros(rank_est,4);
params(:,1)=abs(amps);
params(:,2)=imag(svd_params)/(2*pi);
params(:,3)=real(svd_params)/pi;
params(:,4)=angle(amps);
