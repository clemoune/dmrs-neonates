%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasage des spectres individuels
% JV september 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spec_water fid_phased fid_raw STD_PHI]=phaser_wat(fid_matrix)

global spec_ref spec_p_lb fid_p_lb time factor_for_phi NAA_inf NAA_sup

nb_pts_cplx=size(fid_matrix,1);

dw=1/3012;
lb=3; % lb for phase correction only, not on the final fid

NAA_inf=900; 
NAA_sup=1150;

NS=size(fid_matrix,2);

%%% Generating the raw sum FID %%%%%%%%%%%%

fid_raw=sum(fid_matrix,2);

%%% Generating the reference spectrum for correction %%%%%%%%%%%%

time=((0:nb_pts_cplx-1)*dw)';

fid_ref=fid_raw.*exp(-lb*time)/NS;

spec_ref=fftshift(fft(fid_ref));
spec_ref=spec_ref(NAA_inf:NAA_sup);
spec_water=zeros(NS,length(spec_ref));


fid_phased=zeros(nb_pts_cplx,1);
list_phase=zeros(NS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual spectrum correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options=optimset;
optjv=optimset(options,'MaxFunEvals',1e6);

for p=1:NS
        
        fid_p=squeeze(fid_matrix(:,p));
        
    	fid_p_lb=fid_p.*exp(-lb*time);
	
    	%%% Frequency correction %%% 
            
        y=lsqnonlin('cout_nu',[0 1],[-80 0],[80 3],optjv);
        delta_nu=y(1);
        factor_for_phi=y(2);
	
        correction=exp(2*pi*1i*delta_nu*time);
        
        fid_p=fid_p.*correction;

    	%%% Phase correction %%% 
        
        fid_p_lb=fid_p.*exp(-lb*time);	
        spec_p_lb=fftshift(fft(fid_p_lb));
        spec_p_lb=spec_p_lb(NAA_inf:NAA_sup);
	
        phi=lsqnonlin('cout_phase_spec',0,-pi,pi,optjv);
    	
        list_phase(p)=phi*180/pi;
	
        fid_p=exp(1i*phi)*fid_p;
        
        fid_phased=fid_phased + fid_p;
        
        spec_water(p,:)=spec_ref;

end



MEAN_PHI=mean(list_phase)
STD_PHI=std(list_phase)


fid_raw=fid_raw';
fid_phased=fid_phased';
