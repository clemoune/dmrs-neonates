%%
clearvars -except FolderNames Data_Dir Number_Acq protocol Name_Acquisition data_nb_dir fold
%% Definitions


directory= Data_Dir;
nb_points=2048;
nb_av=16;
nb_coils=4;
Acquisitions=[];
Acquisitions_water =[];
b_values=[];
Acquisitions_water =[];

ProcessedFolder=strcat(directory,filesep,'Processed',filesep);
mkdir(ProcessedFolder)

Number_Acq=size(protocol,2);

for i=1:Number_Acq
    if ischar(protocol(i).type)==1
        
    if ((matches(protocol(i).type,'User:cl_STELASER_PA360_b')==1) && (protocol(i).TM==50) && (protocol(i).Nb_Repetitions==18) && matches(protocol(i).WaterSup,'VAPOR')==1)      
 Acquisitions=[Acquisitions i]; 
 b_values=double(string(split(extractBefore(protocol(i).b_value,' @'))));
 Nb_Repetitions_per_b=double(string(split(extractBefore(protocol(i).rep_per_b,' @'))));
    end

        
        if ((matches(protocol(i).type,'User:cl_STELASER_PA360_b')==1) && (protocol(i).TM==50) && (protocol(i).Nb_Repetitions==18) && matches(protocol(i).WaterSup,'NO_SUPPRESSION')==1)      
        Acquisitions_water=[Acquisitions_water i]; 
        end        
    end
end


%% Uploading files and phasing from different coils


%%% Uploading water


    fid_water=strcat(directory,filesep,num2str(Acquisitions_water),filesep,'fid');    
if exist(eval('fid_water')) == 2
    fid_water=fid_water;
else
    fid_water=strcat(directory,filesep,num2str(Acquisitions_water),filesep,'ser');
end  

[FID_water]=load_array_FID2(fid_water,protocol(Acquisitions_water).Nb_Repetitions);

figure 
    for i=1:size(b_values,1)

   [~,water_ph,~,~]=phaser_wat(FID_water(:,sum(Nb_Repetitions_per_b(1:i-1))+1:sum(Nb_Repetitions_per_b(1:i))));
   spec_ind_coil=abs(fftshift(fft(water_ph)));
   [~,ind_coil]=max(spec_ind_coil);
    
    water_ph1=water_ph.*exp(-1i*2*pi*(ind_coil-1024)*[1:2048]/2048);
        
   Serie_water(i,:) = water_ph1; 
   plot(abs(fftshift(fft(Serie_water(i,:)))))
   hold on
    end
hold off

%%% Establishing the phase list for the coil combination

fid_reading=strcat(directory,filesep,num2str(Acquisitions_water),filesep,'rawdata.job0'); 

if exist(eval('fid_reading')) == 2
[FID1w FID2w FID3w FID4w]=coil_decomposition(fid_reading,nb_points,sum(Nb_Repetitions_per_b),nb_coils);
[~, ~, ~, ~, list_phw]=phaser_wat_recomb([sum(FID1w(:,1:2),2) sum(FID2w(:,1:2),2) sum(FID3w(:,1:2),2) sum(FID4w(:,1:2),2)]);
end 

