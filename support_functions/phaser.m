%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasage des spectres individuels
% JV september 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spec_metab fid_phased fid_raw STD_PHI]=phaser(fid_matrix)

global spec_ref spec_p_lb fid_p_lb time factor_for_phi NAA_inf NAA_sup

nb_pts_cplx=size(fid_matrix,1);

dw=1/3012;
lb=5; % lb for phase correction only, not on the final fid

% Range for upfield metabolites (7T, BW=3012Hz, 2048 pts)

NAA_inf=400;
NAA_sup=800;
NS=size(fid_matrix,2);

%%% Generating the raw sum FID %%%%%%%%%%%%

fid_raw=sum(fid_matrix,2);

%%% Generating the reference spectrum for correction %%%%%%%%%%%%

time=((0:nb_pts_cplx-1)*dw)';

fid_ref=fid_raw.*exp(-lb*time)/NS;
% fid_ref=fid_matrix(:,2).*exp(-lb*time);
spec_ref=fftshift(fft(fid_ref));
spec_ref=spec_ref(NAA_inf:NAA_sup);
spec_metab=zeros(NS,nb_pts_cplx);

fid_phased=zeros(nb_pts_cplx,1);
list_phase=zeros(NS,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Individual spectrum correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options=optimset;
optjv=optimset(options,'MaxFunEvals',1e6,'display', 'off');

for p=1:NS
        
        fid_p=squeeze(fid_matrix(:,p));
        
    	fid_p_lb=fid_p.*exp(-lb*time);
	
    	%%% Frequency correction %%% 
        warning off all    
        y=lsqnonlin('cout_nu',[0 1],[-80 0],[80 3],optjv);
        warning on all   
        delta_nu=y(1);
        factor_for_phi=y(2);
	
        correction=exp(2*pi*1i*delta_nu*time);
        
        fid_p=fid_p.*correction;

    	%%% Phase correction %%% 
        
        fid_p_lb=fid_p.*exp(-lb*time);	
        spec_p_lb=fftshift(fft(fid_p_lb));
        spec_p_lb=spec_p_lb(NAA_inf:NAA_sup);
        
        warning off all 
        phi=lsqnonlin('cout_phase_spec',0,-pi,pi,optjv);
    	warning on all 
        
        list_phase(p)=phi*180/pi;
	
        fid_p=exp(1i*phi)*fid_p;
        
        %%% Frequency correction %%% 
        warning off all     
        y=lsqnonlin('cout_nu',[0 1],[-80 0],[80 3],optjv);
        warning on all 
        delta_nu=y(1);
        factor_for_phi=y(2);
	
        correction=exp(2*pi*1i*delta_nu*time);
        
        fid_p=fid_p.*correction;

        fid_phased=fid_phased + fid_p;

        spec_metab(p,:)=fid_p;

end


MEAN_PHI=mean(list_phase);
STD_PHI=std(list_phase);


fid_raw=fid_raw';
fid_phased=fid_phased';
