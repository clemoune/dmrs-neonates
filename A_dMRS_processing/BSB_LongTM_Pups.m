%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script processes neonates data at long mixing times. 
% The script 'BSB_PhaseEstimated.m' should be run before.
% 1 - Definitions: data preparation
% 2 - Uploading water
% 3 - Uploading metabolites and sorting out outliers (MAD(bv))
% 4 - Save to LCModel without EC correction
% 5 - Correct for EC and/or water residual when possible and ssaving to LCModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clearvars -except Data_Dir Number_Acq protocol Name_Acquisition list_phw fold FolderNames data_nb_dir
%% 1 - Definitions

directory= Data_Dir;
nb_points=2048;
nb_coils=4;
Acquisitions=[];
% b_values=zeros;
MixTime=[];
Acquisitions_water =[];

Nb_Averages=[];


ProcessedFolder=strcat(directory,filesep,'Processed',filesep);
load(strcat(ProcessedFolder,filesep,'ReorderedData.mat'))

Number_Acq=size(protocol,2);

for i=1:Number_Acq
    if ischar(protocol(i).type)==1
        
    if ((matches(protocol(i).type,'User:cl_STELASER_PA360_b')==1) && (protocol(i).TM>100) && (contains(protocol(i).Macro,'Off')==1) && (protocol(i).Nb_Repetitions>6) && matches(protocol(i).WaterSup,'VAPOR')==1)      
 Acquisitions=[Acquisitions i]; 
 b_values(i).b=double(string(split(extractBefore(protocol(i).b_value,' @'))));

 Nb_Repetitions_per_b(i).rep=double(string(split(extractBefore(protocol(i).rep_per_b,' @'))));
 if isnan(Nb_Repetitions_per_b(i).rep)==1
     pre_valTM=split(extractBefore(extractAfter(protocol(i).rep_per_b,'@'),' @'),'*');
    Nb_Repetitions_per_b(i).rep=double(sscanf(string(pre_valTM(2)),'(%d)'))*ones(double(string(pre_valTM(1))),1)
 end
 
 MixTime=[MixTime protocol(i).TM];
    end


        if ((matches(protocol(i).type,'User:cl_STELASER_PA360_b')==1) && (protocol(i).TM>100) && matches(protocol(i).WaterSup,'NO_SUPPRESSION')==1)      
        Acquisitions_water=[Acquisitions_water i]; 
         b_values(i).b=double(string(split(extractBefore(protocol(i).b_value,' @'))));

 Nb_Repetitions_per_b(i).rep=double(string(split(extractBefore(protocol(i).rep_per_b,' @'))));
 if isnan(Nb_Repetitions_per_b(i).rep)==1
     pre_valTM=split(extractBefore(extractAfter(protocol(i).rep_per_b,'@'),' @'),'*');
    Nb_Repetitions_per_b(i).rep=double(sscanf(string(pre_valTM(2)),'(%d)'))*ones(double(string(pre_valTM(1))),1)
 end
 
        end 
    end
end



%% 2 - Uploading water

for k=1:size(Acquisitions_water,2)

    fid_water=strcat(directory,filesep,num2str(Acquisitions_water(1,k)),filesep,'fid');    
if exist(eval('fid_water')) == 2
    fid_water=fid_water;
else
    fid_water=strcat(directory,filesep,num2str(Acquisitions_water(1,k)),filesep,'ser');
end  

[FID_water]=load_array_FID2(fid_water,sum(Nb_Repetitions_per_b(Acquisitions_water(1,k)).rep));

    for i=1:size(b_values(Acquisitions_water(1,k)).b,1)

   [spec_water water_ph water_raw STD_PHI]=phaser_wat(FID_water(:,sum(Nb_Repetitions_per_b(Acquisitions_water(1,k)).rep(1:i-1))+1:sum(Nb_Repetitions_per_b(Acquisitions_water(1,k)).rep(1:i))));
    spec_ind_coil=abs(fftshift(fft(water_ph)));

    % centering the water peak 
    [M,ind_coil]=max(spec_ind_coil);
    water_ph1=water_ph.*exp(-1i*2*pi*(ind_coil-1024)*[1:2048]/2048);
        
    WaterTM(Acquisitions_water(1,k)).Serie_water(i,:) = water_ph1;  
    end   

    
    b0=find(b_values(Acquisitions(k)).b(1:end)<3000);
    b3=find(b_values(Acquisitions(k)).b(1:end)>3000);
    [spec_water b0_spec b0_spec_raw STD_PHI]=phaser(WaterTM(Acquisitions_water(1,k)).Serie_water(b0,:)');
    [spec_water b3_spec b3_spec_raw STD_PHI]=phaser(WaterTM(Acquisitions_water(1,k)).Serie_water(b3,:)');
    WaterTM(Acquisitions_water(1,k)).Serie_water_tot(1,1:2048)=b0_spec;
    WaterTM(Acquisitions_water(1,k)).Serie_water_tot(2,1:2048)=b3_spec;
    
    figure 
    plot(abs(fftshift(fft(WaterTM(Acquisitions_water(1,k)).Serie_water_tot(1,1:2048)))))
    hold on 
    plot(abs(fftshift(fft(WaterTM(Acquisitions_water(1,k)).Serie_water_tot(2,1:2048)))))
    hold off
end


name_water=strcat(ProcessedFolder,filesep,'Water_TM','.mat');
save(name_water,'WaterTM','-mat')



%% 3 - Uploading metabolites and sorting out outliers (MAD(bv))
% only b=0 and a low b-value (similar MAD at 2 b-values)
% we also do not make a distinction as function of direction 

for i=1:3%size(Acquisitions,2)
    
fid_reading=strcat(directory,filesep,num2str(Acquisitions(i)),filesep,'rawdata.job0'); 

if exist(eval('fid_reading')) == 2

[FID1 FID2 FID3 FID4]=coil_decomposition(fid_reading,nb_points,16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep),nb_coils);
Recombined_FID=FID1*exp(1i*list_phw(1)*pi/180)+FID2*exp(1i*list_phw(2)*pi/180)+FID3*exp(1i*list_phw(3)*pi/180)+FID4*exp(1i*list_phw(4)*pi/180);



lb =10;

range=[390:940];
range_downfield=[390:940]+1000;

% range=[600:800];
% range_downfield=[600:800]+1000;

mean_signal=[];
mean_downfield=[];
median_spec=[];
median_plot=[];
std_downfield=[];
    
std_median_spec=[];
std_median_plot=[];

% calculating Mup(i), Mdownfield and STDdownfield

for j=1:size(Recombined_FID,2)
   spec_abs= abs(fftshift(fft(Recombined_FID(:,j).*exp(-lb*[1:2048]'/3012))));
   mean_signal=[mean_signal sum(spec_abs(range))/size(range,2)];
   mean_downfield=[mean_downfield sum(spec_abs(range_downfield))/size(range,2)];
   std_downfield=[std_downfield std(spec_abs(range_downfield))];
   
end

for k=1:size(Nb_Repetitions_per_b(Acquisitions(i)).rep,1)
    median_spec=[median_spec median(mean_signal(1,16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:k-1))+1:16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:k))))];
    std_median_spec=[std_median_spec std(mean_signal(1,16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:k-1))+1:16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:k))))/mean(mean_signal(1,16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:k-1))+1:16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:k))))];
end

for kk=1:size(Nb_Repetitions_per_b(Acquisitions(i)).rep,1)
median_plot=[median_plot;ones(Nb_Repetitions_per_b(Acquisitions(i)).rep(k)*16,1)*median_spec(kk)];
std_median_plot=[std_median_plot;ones(Nb_Repetitions_per_b(Acquisitions(i)).rep(k)*16,1)*std_median_spec(kk)];
end

[keepers]= find((mean_signal'-median_plot.*(1-std_median_plot))>0)

Serie_metab(i).keepers=keepers;
Serie_metab(i).median_spec=median_spec;
Serie_metab(i).std_median_spec=std_median_spec;
Serie_metab(i).SNR=mean_signal./std_downfield;

figure
plot(mean_signal)
hold on 
plot(mean_downfield)
plot(median_plot, 'k', 'LineWidth', 1)
plot(median_plot.*(1-std_median_plot), 'r', 'LineWidth', 1)
plot(median_plot.*(1+std_median_plot), 'r', 'LineWidth', 1)


    for j=1:size(Nb_Repetitions_per_b(Acquisitions(i)).rep,1)

            
        BRepetKept=find(16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:j-1))+1 <= keepers & keepers <= 16*sum(Nb_Repetitions_per_b(Acquisitions(i)).rep(1:j)));
        Serie_metab(i).BRepetKept(j).reps=Serie_metab(i).keepers(BRepetKept);
                
        [~,fid_phasedtot,fid_rawtot,~]=phaser(Recombined_FID(:,Serie_metab(i).keepers(BRepetKept)));


%      fid_phasedtot=fid_phasedtot.*exp(-1i*2*pi*(982-1024)*[1:2048]/2048)

%     spec_ind_coil=abs(fftshift(fft(fid_phasedtot)));
%     [M,ind_coil]=max(spec_ind_coil(1,1250:1350));
%     ind_coil=ind_coil+1250;    
%     fid_phasedtot=fid_phasedtot.*exp(-1i*2*pi*(ind_coil-1330)*[1:2048]/2048);
    

    Serie_metab(i).spec(j,:)=fid_phasedtot/(size(BRepetKept,1)*nb_coils);

       end       

end
end

for i=1:3%size(Acquisitions,2)    
    b0=find(b_values(Acquisitions(i)).b(1:end)<3000);
    b3=find(b_values(Acquisitions(i)).b(1:end)>3000);
    [spec_water b0_spec b0_spec_raw STD_PHI]=phaser(Serie_metab(i).spec(b0,:)');
    [spec_water b3_spec b3_spec_raw STD_PHI]=phaser(Serie_metab(i).spec(b3,:)');
    Serie_metab_tot(i).spec(1,1:2048)=b0_spec;
    Serie_metab_tot(i).spec(2,1:2048)=b3_spec;
    figure 
    plot(abs(fftshift(fft(Serie_metab_tot(i).spec(1,1:2048)))))
    hold on 
    plot(abs(fftshift(fft(Serie_metab_tot(i).spec(2,1:2048)))))
    hold off
end



name_serie_acq=strcat(ProcessedFolder,filesep,'LongTMData','.mat');
save(name_serie_acq,'Serie_metab_tot','Serie_metab','-mat')

%% 4 - Save to LCModel without EC correction, or residual water removal

%  no eddy current correction 

for i = 1:3%size(Acquisitions,2)
for k=1:2  

directory_LCModel_up=strcat('Fitting_LCModel_202210');  

directory_LCModelBV = strcat(directory_LCModel_up,filesep,Name_Acquisition,filesep,'TM',num2str(MixTime(i)),filesep,'b',num2str(k),filesep,'BV',filesep);
mkdir(directory_LCModelBV)
save_fid(Serie_metab_tot(i).spec(k,:)',0,directory_LCModelBV);
   
end

end

%% 5 - Correct for EC and/or water residual when possible and ssaving to LCModel
% one spectrum per b value - eddy current correction - 
try
for k=1:3
    for bi=1:2
        
%%% Performing EC correction & water residual removal

%     [td_synth, ~, ~]=svdfid(WaterTM(Acquisitions_water(k)).Serie_water_tot(bi,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
%     fid_cor=deconv_ECC(Serie_metab_tot(k).spec(bi,:)',td_synth);
%     [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);


%%% Performing water residual removal only

[~, td_diff, ~]=svdfid(Serie_metab_tot(k).spec(bi,:)', 8, 3012, -150, 150, 3, 0.001, 2000);
 



        Serie_metab_tot(k).spec_EC(bi,:)=td_diff';
%       Serie_metab_tot(k).spec_EC(bi,:)=Serie_metab_tot(k).spec(bi,:);

    figure
    hold on 
    plot(real(fftshift(fft(Serie_metab_tot(k).spec(bi,:)))))
    plot(real(fftshift(fft(td_diff'))))
    hold off
   
end

end

catch
for k=1:3
    for bi=1
        
%%% Performing EC correction & water residual removal

%     [td_synth, ~, ~]=svdfid(WaterTM(Acquisitions_water(k)).Serie_water_tot(bi,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
%     fid_cor=deconv_ECC(Serie_metab_tot(k).spec(bi,:)',td_synth);
%     [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);


%%% Performing water residual removal only

        [~, td_diff, ~]=svdfid(Serie_metab_tot(k).spec(bi,:)', 8, 3012, -150, 150, 3, 0.001, 2000);
 

        
        Serie_metab_tot(k).spec_EC(bi,:)=td_diff';
%       Serie_metab_tot(k).spec_EC(bi,:)=Serie_metab_tot(k).spec(bi,:);

    figure
    hold on 
    plot(real(fftshift(fft(Serie_metab_tot(k).spec(bi,:)))))
    plot(real(fftshift(fft(td_diff'))))
    hold off
   
end

end

% for k=2
%     for bi=1

for k=1:3
    for bi=2

%%% Performing EC correction & water residual removal

%     [td_synth, ~, ~]=svdfid(WaterTM(Acquisitions_water(k)).Serie_water_tot(bi,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
%     fid_cor=deconv_ECC(Serie_metab_tot(k).spec(bi,:)',td_synth);
%     [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);

%%% Performing water residual removal only

%        [~, td_diff, ~]=svdfid(Serie_metab_tot(k).spec(bi,:)', 8, 3012, -150, 150, 3, 0.001, 2000);
 

%         Serie_metab_tot(k).spec_EC(bi,:)=td_diff';

%%% Neither performing EC corr, nor water residual suppression

      Serie_metab_tot(k).spec_EC(bi,:)=Serie_metab_tot(k).spec(bi,:);

    figure
    plot(real(fftshift(fft(Serie_metab_tot(k).spec(bi,:)))))

   
    end
end
end

name_serie_acq=strcat(ProcessedFolder,filesep,'LongTMData','.mat');
save(name_serie_acq,'Serie_metab_tot','Serie_metab','-mat')

% Save the corrected spectra to LCModel 
for i = 1:3%size(Acquisitions,2)
for k=1:2  

directory_LCModelBV = strcat(directory_LCModel_up,filesep,'EC_',Name_Acquisition,filesep,'TM',num2str(MixTime(i)),filesep,'b',num2str(k),filesep,'BV',filesep);
mkdir(directory_LCModelBV)
save_fid(Serie_metab_tot(i).spec_EC(k,:)',0,directory_LCModelBV);  
   
end

end
