%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script processes neonates data at high-b values. 
% It notably generates a big structure called "Reordered" saving all the
% processing steps. It can probably be greatly simplified / generalised.
% Script structure:
% 1 - Definitions: preparing data.
% 2 - Uploading files, combining coils.
% 3 - Calculating the mean upfield Mup(i) signal for each transient i
% 4 - Estimating diffusion & kurtosis tensors from all data
% 5 - Calculating the median and MAD at the b-value level
% 6A - Outlier detection | median & MAD(bv)
%   this first step of outliers detection is to remove gross outliers
%   and re-evaluate the kurtosis tensor without those (to avoid corruption 
%   by motion with a preferential direction)
% 6A-bis Computing the FA of the mean upfield signal to check for 
% directionality dependence
% 6B - Outlier detection | median & MAD(bv=0). For isotropic voxels (FA~0).
%   Regardless of the gradient direction, the MAD calculated at b=0 is 
%   applied around the median of Mup at each bval to detect outliers. 
%   The MAD should be the same at all b-values: an increase in MAD with bv
%   underlines contribution of motion. However this approach only works in
%   isotropic voxels: as soon as there is a direction dependency, the
%   DKI tensor estimate should be used. 
% 6C - Outlier detection based on the DKI prediction | median & MAD(bv=0)
% 6D - Outlier detection based on the DKI prediction following 6A | median & MAD(bv=0) 
% 7 - Summing spectra
% 8 - Save to LCModel prior to eddy current correction & water residual removal
% 9 - Eddy current correction & save to LCModel post EC corr & water residual removal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clearvars -except FolderNames Data_Dir Number_Acq protocol Name_Acquisition data_nb_dir fold

%% 1 - Definitions

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


%% 2 - Uploading files, combining coils.


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

name_water=strcat(ProcessedFolder,filesep,'water_high_b','.mat');
save(name_water,'Serie_water','list_phw','-mat')

%%% Uploading metabolites and applying coil combination results

% There are two acquisitions of metabolites (to be averaged)

for j=1:size(Acquisitions,2)
    
fid_reading=strcat(directory,filesep,num2str(Acquisitions(j)),filesep,'rawdata.job0'); 

if exist(eval('fid_reading')) == 2

[FID1 FID2 FID3 FID4]=coil_decomposition(fid_reading,nb_points,nb_av*sum(Nb_Repetitions_per_b),nb_coils);
AllRaw(j).Recombined_FID=FID1*exp(1i*list_phw(1)*pi/180)+FID2*exp(1i*list_phw(2)*pi/180)+FID3*exp(1i*list_phw(3)*pi/180)+FID4*exp(1i*list_phw(4)*pi/180);

end
end 

% Sometimes, there is a frequency difference between the 2 rounds that is
% too important for the phaser (usually due to a temperature drift). We
% look for this frequency difference and adjust.

Spec1=sum(AllRaw(1).Recombined_FID(:,1:16),2);
Spec2=sum(AllRaw(2).Recombined_FID(:,1:16),2);
Spec1=Spec1';
Spec2=Spec2';
Spec1Abs=abs(fftshift(fft(Spec1)));
Spec2Abs=abs(fftshift(fft(Spec2)));
NAA_range=[1450:1700];

[m ind1]=max(Spec1Abs(NAA_range));
[m ind2]=max(Spec2Abs(NAA_range));

if ind1 > ind2   
AllRaw(2).Recombined_FID=(AllRaw(2).Recombined_FID).*exp(1i*(ind2-ind1)*2*pi*[1:2048]'/2048);
elseif ind2 > ind1
AllRaw(1).Recombined_FID=(AllRaw(1).Recombined_FID).*exp(1i*(ind1-ind2)*2*pi*[1:2048]'/2048);    
end

%%% Sorting out by b-value and direction

TotalRecombined=[AllRaw(1).Recombined_FID AllRaw(2).Recombined_FID];

for bv=1:6
    for dir=1:16
        Reordered(bv).Direction(dir).indices=[];
        for i=1:Nb_Repetitions_per_b(bv)
        Reordered(bv).Direction(dir).indice(i).list= [sum(Nb_Repetitions_per_b(1:bv-1))*nb_av+dir+(i-1)*nb_av sum(Nb_Repetitions_per_b(1:bv-1))*nb_av+dir+(i-1)*nb_av+288];
        Reordered(bv).Direction(dir).indices=[Reordered(bv).Direction(dir).indices Reordered(bv).Direction(dir).indice(i).list];
        end
         Reordered(bv).Direction(dir).spectra=TotalRecombined(:,Reordered(bv).Direction(dir).indices);
    end
end


clearvars -except FolderNames data_nb_dir fold Data_Dir Number_Acq protocol Name_Acquisition Reordered Serie_water Acquisitions Acquisitions_water b_values nb_av nb_coils nb_points ProcessedFolder directory Nb_Repetitions_per_b


%% 3 - Calculating the mean upfield Mup(i) signal for each transient i
% also calculating the downfield signal (approximation of noise) and its STD 
% derive SNR/transient


range=[390:940];
range_downfield=[390:940]+1000;

% initiating variables
for bv=1:6
    for dir=1:nb_av
           Reordered(bv).Direction(dir).mean_signal=[];
           Reordered(bv).Direction(dir).mean_downfield=[];
           Reordered(bv).Direction(dir).std_downfield=[];
    end
    
Reordered(bv).Spectra=[ ];
Reordered(bv).MeanSignal=[ ];
Reordered(bv).MeanSignalDownfield=[];

end


ReorderedMean=[ ];
ReorderedMeanDownfield=[ ];


%%% Calculating the mean signal and storing values at the direction &
%%% b-value levels

for bv=1:6
    for dir=1:nb_av
    for i=1:size(Reordered(bv).Direction(dir).spectra,2)  
    spec_abs= abs(fftshift(fft(Reordered(bv).Direction(dir).spectra(:,i).*exp(-10*[1:2048]'/3012))));
    Reordered(bv).Direction(dir).mean_signal=[Reordered(bv).Direction(dir).mean_signal sum(spec_abs(range))/size(range,2)];
    Reordered(bv).Direction(dir).mean_downfield=[Reordered(bv).Direction(dir).mean_downfield sum(spec_abs(range_downfield))/size(range,2)];
    Reordered(bv).Direction(dir).std_downfield=[Reordered(bv).Direction(dir).std_downfield std(spec_abs(range_downfield))];
    end
    Reordered(bv).Direction(dir).SNRrel=Reordered(bv).Direction(dir).mean_signal./Reordered(bv).Direction(dir).std_downfield;
    end
end


 
for bv=1:6
       %concatenating per direction   
    for dir=1:nb_av
    Reordered(bv).Spectra=[Reordered(bv).Spectra Reordered(bv).Direction(dir).spectra];
    Reordered(bv).MeanSignal=[Reordered(bv).MeanSignal Reordered(bv).Direction(dir).mean_signal];
    Reordered(bv).MeanSignalDownfield=[Reordered(bv).MeanSignalDownfield Reordered(bv).Direction(dir).mean_downfield];

    end 
    ReorderedMean=[ReorderedMean Reordered(bv).MeanSignal];
    ReorderedMeanDownfield=[ReorderedMeanDownfield Reordered(bv).MeanSignalDownfield];
end 



%% 4 - Estimating diffusion & kurtosis tensors from all data

%%% Gradient lists extracted from the Bruker method file.
Gdiff_list1(1,:)=[-0.513259950761164 1.03215602988595 1.39814965668857 -0.0540160779674595 -0.732906813907144 -1.69551358163122 -1.91381573739913 1.05151790945474 2.04966481351727 1.37215233237696 -0.465645572446854 -0.54129669265212 -1.61899477395015 1.062752250192 0.384753945731602 -0.814043418560211];
Gdiff_list1(1,:)=Gdiff_list1(1,:)./norm(Gdiff_list1(1,:));
Gdiff_list2(1,:)=[0.0987780515121413 1.04967641703115 -0.650268798983989 -1.60676351131909 1.62374528049372 -0.74468760784915 0.779969377657589 1.77457721662268 0.161302602950086 -1.51426768199829 -1.98312418230019 1.8132360288041 -0.839268226643984 0.643494301139662 -1.03214043844737 0.429524024760702];
Gdiff_list2(1,:)=Gdiff_list2(1,:)./norm(Gdiff_list2(1,:));
Gdiff_list3(1,:)=[1.99972185425825 1.45083208001915 1.37637509476287 1.29902796618415 1.0480358968646 0.918016577265322 -0.0322097518431654 -0.131400773413133 -0.211981125161449 -0.310276874616744 -0.349961502669672 -0.831416483649785 -0.972915793160905 -1.65183325854233 -1.74892278030781 -1.8506546949029];
Gdiff_list3(1,:)=Gdiff_list3(1,:)./norm(Gdiff_list3(1,:));


%%% Using all data to fit the tensor

Gdiff1=[];
Gdiff2=[];
Gdiff3=[];

for bv=1:6
for dir=1:16
    Gdiff1=[Gdiff1 repmat(Gdiff_list1(1,dir),1,sum(size(Acquisitions,2)*Nb_Repetitions_per_b(bv)))];
    Gdiff2=[Gdiff2 repmat(Gdiff_list2(1,dir),1,sum(size(Acquisitions,2)*Nb_Repetitions_per_b(bv)))];
    Gdiff3=[Gdiff3 repmat(Gdiff_list3(1,dir),1,sum(size(Acquisitions,2)*Nb_Repetitions_per_b(bv)))];    
end
end

bvecs=[Gdiff1;Gdiff2;Gdiff3];
bvals=[];
Data=[];

for bv=1:6
bvals=[bvals ones(1,nb_av*size(Acquisitions,2)*Nb_Repetitions_per_b(bv))*b_values(bv)];
Data=[Data Reordered(bv).MeanSignal];
end

% Calculating DKI and DTI for all possible b-max
Reordered(2).DTI=diffmodelfit(Data(1,1:128),bvecs(:,1:128),bvals(1,1:128),'DTI');
Reordered(3).DTI=diffmodelfit(Data(1,1:224),bvecs(:,1:224),bvals(1,1:224),'DTI');
Reordered(4).DTI=diffmodelfit(Data(1,1:224+96),bvecs(:,1:224+96),bvals(1,1:224+96),'DTI');
Reordered(5).DTI=diffmodelfit(Data(1,1:224+96+128),bvecs(:,1:224+96+128),bvals(1,1:224+96+128),'DTI');
Reordered(6).DTI=diffmodelfit(Data,bvecs,bvals,'DTI');


Reordered(2).DKI=diffmodelfit(Data(1,1:128),bvecs(:,1:128),bvals(1,1:128),'DKI');
Reordered(3).DKI=diffmodelfit(Data(1,1:224),bvecs(:,1:224),bvals(1,1:224),'DKI');
Reordered(4).DKI=diffmodelfit(Data(1,1:224+96),bvecs(:,1:224+96),bvals(1,1:224+96),'DKI');
Reordered(5).DKI=diffmodelfit(Data(1,1:224+96+128),bvecs(:,1:224+96+128),bvals(1,1:224+96+128),'DKI');
Reordered(6).DKI=diffmodelfit(Data,bvecs,bvals,'DKI');


figure
title('DTI predictions')
for i =2:6
plot(ReorderedMean,'k')
hold on 
plot(Reordered(i).DTI.pred)
end

figure
title('DKI predictions')
for i =2:6
plot(ReorderedMean,'k')
hold on 
plot(Reordered(i).DKI.pred)
end

% We plot the DKI estimation where all b-vals are used, it is ultimately the one kept for the prediction.
DKI_pred=figure
plot(ReorderedMean,'k')
hold on 
plot(Reordered(6).DKI.pred)
name_fig=strcat(ProcessedFolder,filesep,'DKIpred.fig');
savefig(DKI_pred, name_fig)


for bv=1:6
    for dir=1:nb_av
    Reordered(bv).Direction(dir).pred=Reordered(6).DKI.pred(1,sum(Nb_Repetitions_per_b(1:bv-1)*2*nb_av)+ 1 + (dir-1)*Nb_Repetitions_per_b(bv)*2:sum(Nb_Repetitions_per_b(1:bv-1)*2*nb_av)+ (dir)*Nb_Repetitions_per_b(bv)*2);
    end 
end 

%% 5 - Calculating the median and MAD at the b-value level

% arbitrary factor to define the outliers threshold, usually 3xsigma.
FactorMAD=3;

%colours for histograms
coloursHist=[226 240 217; 197 224 180; 169 209 142; 84 130 53; 56 87 35; 24 38 16]/256;

%initiating variables
Median_plot=[];
 
MAD_plot.BV=[]; % median signal of Mup - mean absolute deviation defined bval by bval
MAD_plot.BV_plus=[]; % median signal of Mup + mean absolute deviation defined bval by bval

MAD_plot.B0=[]; % median signal of Mup - mean absolute deviation defined at b=0
MAD_plot.B0_plus=[]; % median signal of Mup + mean absolute deviation defined at b=0

MAD_plot.B0max=[]; % median signal pre-sorted of Mup - mean absolute deviation defined at b=0
MAD_plot.B0max_plus=[]; % median signal pre-sorted of Mup + mean absolute deviation defined at b=0




for bv=1:6
    
    %defining the median and MAD over all directions
    Reordered(bv).median_spec=median(Reordered(bv).MeanSignal);
    Reordered(bv).MAD=median(abs(Reordered(bv).MeanSignal-Reordered(bv).median_spec));
    
    Reordered(bv).OutliersBV.SkewnessInit=skewness(Reordered(bv).MeanSignal);
    
    Median_plot=[Median_plot Reordered(bv).median_spec*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    MAD_plot.BV=[MAD_plot.BV (Reordered(bv).median_spec - FactorMAD*Reordered(bv).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    MAD_plot.BV_plus=[MAD_plot.BV_plus (Reordered(bv).median_spec + FactorMAD*Reordered(bv).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
end 


%% 6A - Outlier detection | median & MAD(bv)
% this first step of outliers detection is to remove gross motion outliers
% and re-evaluate the kurtosis tensor without those (to avoid corruption by motion with a preferential direction)

for bv=1:6
Reordered(bv).WeightedSpectraBV=[ ];
Reordered(bv).OutliersBV.SortedMeanSignal=[ ];
Reordered(bv).OutliersBV.TotalKept=[];
end

ReorderedMeanBV=[];
    bvecsBV1=[];
    bvecsBV2=[];
    bvecsBV3=[];
    bvalsBV=[];
    bvecsBV=[];

for bv=1:6

    % sorting out at the direction & b-value levels (but threshold defined at the b-value level)
    EjectedDirCounter=0;
    for dir=1:nb_av
        % defining the transients kept by directions
        Reordered(bv).Direction(dir).keepers= find((Reordered(bv).Direction(dir).mean_signal-Reordered(bv).median_spec+FactorMAD*Reordered(bv).MAD)>0);

    
    if size(Reordered(bv).Direction(dir).keepers,2)>0
    Reordered(bv).Direction(dir).spectraBV=Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepers);
    Reordered(bv).Direction(dir).mean_signalBV=Reordered(bv).Direction(dir).mean_signal(:,Reordered(bv).Direction(dir).keepers);        
    Reordered(bv).Direction(dir).weightedspectraBV=Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepers)/size(Reordered(bv).Direction(dir).keepers,2);
    Reordered(bv).WeightedSpectraBV=[Reordered(bv).WeightedSpectraBV Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepers)/size(Reordered(bv).Direction(dir).keepers,2)];
    Reordered(bv).OutliersBV.SortedMeanSignal=[Reordered(bv).OutliersBV.SortedMeanSignal Reordered(bv).Direction(dir).mean_signal(Reordered(bv).Direction(dir).keepers)];
    Reordered(bv).Direction(dir).nbkeeper=size(Reordered(bv).Direction(dir).keepers,2);
    ReorderedMeanBV=[ReorderedMeanBV Reordered(bv).Direction(dir).mean_signal(Reordered(bv).Direction(dir).keepers)];
    bvecsBV1=[bvecsBV1 Gdiff_list1(1,dir)*ones(1,size(Reordered(bv).Direction(dir).keepers,2))];
    bvecsBV2=[bvecsBV2 Gdiff_list2(1,dir)*ones(1,size(Reordered(bv).Direction(dir).keepers,2))];
    bvecsBV3=[bvecsBV3 Gdiff_list3(1,dir)*ones(1,size(Reordered(bv).Direction(dir).keepers,2))];
    bvalsBV=[bvalsBV b_values(bv)*ones(1,size(Reordered(bv).Direction(dir).keepers,2))];
    else 
        % eject the whole direction if no transient pass the threshold
    EjectedDirCounter= EjectedDirCounter + 1;
    Reordered(bv).WeightedSpectraBV=Reordered(bv).WeightedSpectraBV;
    Reordered(bv).Direction(dir).nbkeeper=0;
    end
     
    Reordered(bv).OutliersBV.TotalKept=[Reordered(bv).OutliersBV.TotalKept Reordered(bv).Direction(dir).nbkeeper];
    end
    Reordered(bv).WeightedSpectraBV=Reordered(bv).WeightedSpectraBV*nb_av/(nb_av-EjectedDirCounter);
    % computing new skewness and % of outliers
    Reordered(bv).OutliersBV.SkewnessSorted=skewness(Reordered(bv).OutliersBV.SortedMeanSignal);
    Reordered(bv).OutliersBV.Percent=1-sum(Reordered(bv).OutliersBV.TotalKept)/size(Reordered(bv).Spectra,2);
    bvecsBV=[bvecsBV1;bvecsBV2;bvecsBV3];


end

%% 6A-bis Computing the FA of the mean upfield signal to check for directionality dependence
% Using the previous estimation where larger outliers have been removed. 

Reordered(2).DKIBV=diffmodelfit(ReorderedMeanBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)),bvecsBV(:,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)),bvalsBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)),'DKI');
Reordered(3).DKIBV=diffmodelfit(ReorderedMeanBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)),bvecsBV(:,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)),bvalsBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)),'DKI');
Reordered(4).DKIBV=diffmodelfit(ReorderedMeanBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)+size(Reordered(4).WeightedSpectraBV,2)),bvecsBV(:,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)+size(Reordered(4).WeightedSpectraBV,2)),bvalsBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)+size(Reordered(4).WeightedSpectraBV,2)),'DKI');
Reordered(5).DKIBV=diffmodelfit(ReorderedMeanBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)+size(Reordered(4).WeightedSpectraBV,2)+size(Reordered(5).WeightedSpectraBV,2)),bvecsBV(:,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)+size(Reordered(4).WeightedSpectraBV,2)+size(Reordered(5).WeightedSpectraBV,2)),bvalsBV(1,1:size(Reordered(1).WeightedSpectraBV,2)+size(Reordered(2).WeightedSpectraBV,2)+size(Reordered(3).WeightedSpectraBV,2)+size(Reordered(4).WeightedSpectraBV,2)+size(Reordered(5).WeightedSpectraBV,2)),'DKI');
Reordered(6).DKIBV=diffmodelfit(ReorderedMeanBV,bvecsBV,bvalsBV,'DKI');


figure
title('DKI predictions')
for i =2:6
plot(ReorderedMeanBV,'k')
hold on 
plot(Reordered(i).DKIBV.pred)
end

DKIBV_pred=figure
plot(ReorderedMeanBV,'k')
hold on 
plot(Reordered(6).DKIBV.pred)
name_fig=strcat(ProcessedFolder,filesep,'DKIBVpred.fig');
savefig(DKIBV_pred, name_fig)

for i=2:6
 Reordered(i).DKIBV.FA
end

counter_predbv=[];
for bv=1:6
    counter_predbv=[counter_predbv sum(Reordered(bv).OutliersBV.TotalKept)];
    counter_pred=[];
    for dir=1:nb_av
        counter_pred=[counter_pred 1+sum(counter_predbv(1:bv-1))+sum(Reordered(bv).OutliersBV.TotalKept(1:dir-1))];
    end 
    Reordered(bv).predBV=Reordered(6).DKIBV.pred(1,counter_pred);
end



%% 6B - Outlier detection | median & MAD(bv=0). For isotropic voxels (FA~0).

for bv=1:6
Reordered(bv).WeightedSpectraB0=[ ];
Reordered(bv).OutliersB0.SortedMeanSignal=[ ];
Reordered(bv).OutliersB0.TotalKept=[];
end

for bv=1:6

    % sorting out at the direction & b-value levels (but threshold defined at the b-value level)

    MAD_plot.B0=[MAD_plot.B0 (Reordered(bv).median_spec - FactorMAD*Reordered(1).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    MAD_plot.B0_plus=[MAD_plot.B0_plus (Reordered(bv).median_spec + FactorMAD*Reordered(1).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];

    EjectedDirCounter=0;  
    
    for dir=1:nb_av
    Reordered(bv).Direction(dir).keepersB0= find((Reordered(bv).Direction(dir).mean_signal-Reordered(bv).median_spec+FactorMAD*Reordered(1).MAD)>0);
    

    
    if size(Reordered(bv).Direction(dir).keepersB0,2)>0
    Reordered(bv).Direction(dir).weightedspectraB0=Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepersB0)/size(Reordered(bv).Direction(dir).keepersB0,2);
    Reordered(bv).WeightedSpectraB0=[Reordered(bv).WeightedSpectraB0 Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepersB0)/size(Reordered(bv).Direction(dir).keepersB0,2)];
    Reordered(bv).OutliersB0.SortedMeanSignal=[Reordered(bv).OutliersB0.SortedMeanSignal Reordered(bv).Direction(dir).mean_signal(Reordered(bv).Direction(dir).keepersB0)];
    Reordered(bv).Direction(dir).nbkeeperB0=size(Reordered(bv).Direction(dir).keepersB0,2);
    else 
    EjectedDirCounter=EjectedDirCounter+1;
    Reordered(bv).WeightedSpectraB0=Reordered(bv).WeightedSpectraB0;
    Reordered(bv).Direction(dir).nbkeeperB0=0;
    end
    Reordered(bv).OutliersB0.TotalKept=[Reordered(bv).OutliersB0.TotalKept Reordered(bv).Direction(dir).nbkeeperB0];    
    end
    
    Reordered(bv).WeightedSpectraB0=Reordered(bv).WeightedSpectraB0*nb_av/(nb_av-EjectedDirCounter);
    % computing new skewness and % of outliers
    Reordered(bv).OutliersB0.SkewnessSorted=skewness(Reordered(bv).OutliersB0.SortedMeanSignal);
    Reordered(bv).OutliersB0.Percent=1-sum(Reordered(bv).OutliersB0.TotalKept)/size(Reordered(bv).Spectra,2);
end


%% 6C - Outlier detection based on the DKI prediction | median & MAD(bv=0)


MAD_plot.DKI=[];
MAD_plot.DKI_plus=[];
DKI_pred_plot=[];

for bv=1:6
Reordered(bv).WeightedSpectraDKI=[ ];
Reordered(bv).OutliersDKI.SortedMeanSignal=[ ];
Reordered(bv).OutliersDKI.TotalKept=[];

Reordered(bv).Pred=[];
end

for bv=1:6
    for dir=1:nb_av
    Reordered(bv).Pred=[Reordered(bv).Pred Reordered(bv).Direction(dir).pred];
    end 
end 


for bv=1
    MAD_plot.DKI=[MAD_plot.DKI (Reordered(bv).median_spec - FactorMAD*Reordered(1).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    MAD_plot.DKI_plus=[MAD_plot.DKI_plus (Reordered(bv).median_spec + FactorMAD*Reordered(1).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    DKI_pred_plot=[DKI_pred_plot Reordered(bv).median_spec*ones(1,nb_av*Nb_Repetitions_per_b(bv)*size(Acquisitions,2))];
    for dir=1:nb_av
    Reordered(bv).Direction(dir).weightedspectraDKI=Reordered(bv).Direction(dir).weightedspectraB0;
    Reordered(bv).WeightedSpectraDKI=Reordered(bv).WeightedSpectraB0;
    
    Reordered(bv).Direction(dir).nbkeeperDKI=Reordered(bv).Direction(dir).nbkeeperB0;
    end
    Reordered(bv).OutliersDKI.SortedMeanSignal=Reordered(bv).OutliersB0.SortedMeanSignal;
    Reordered(bv).OutliersDKI.SkewnessSorted=Reordered(bv).OutliersB0.SkewnessSorted;
    Reordered(bv).OutliersDKI.Percent=Reordered(bv).OutliersB0.Percent;
  
end


for bv=2:6

    % sorting out at the direction & b-value levels (but threshold defined at the b-value level)

     EjectedDirCounter=0;    
     
    for dir=1:nb_av
        
    MAD_plot.DKI=[MAD_plot.DKI (Reordered(bv).Direction(dir).pred-mean(Reordered(bv).Pred)+Reordered(bv).median_spec - FactorMAD*Reordered(1).MAD)];
    MAD_plot.DKI_plus=[MAD_plot.DKI_plus (Reordered(bv).Direction(dir).pred-mean(Reordered(bv).Pred)+Reordered(bv).median_spec + FactorMAD*Reordered(1).MAD)];
    DKI_pred_plot=[DKI_pred_plot Reordered(bv).Direction(dir).pred-mean(Reordered(bv).Pred)+Reordered(bv).median_spec];
    
    Reordered(bv).Direction(dir).keepersDKI= find((Reordered(bv).Direction(dir).mean_signal -(Reordered(bv).Pred(dir)-mean(Reordered(bv).Pred)+Reordered(bv).median_spec)+FactorMAD*Reordered(1).MAD)>0);

    
    if size(Reordered(bv).Direction(dir).keepersDKI,2)>0
    Reordered(bv).Direction(dir).weightedspectraDKI=Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepersDKI)/size(Reordered(bv).Direction(dir).keepersDKI,2);
    Reordered(bv).WeightedSpectraDKI=[Reordered(bv).WeightedSpectraDKI Reordered(bv).Direction(dir).spectra(:,Reordered(bv).Direction(dir).keepersDKI)/size(Reordered(bv).Direction(dir).keepersDKI,2)];
    Reordered(bv).OutliersDKI.SortedMeanSignal=[Reordered(bv).OutliersDKI.SortedMeanSignal Reordered(bv).Direction(dir).mean_signal(Reordered(bv).Direction(dir).keepersDKI)];
    Reordered(bv).Direction(dir).nbkeeperDKI=size(Reordered(bv).Direction(dir).keepersDKI,2);
    else 
    EjectedDirCounter=EjectedDirCounter+1;
    Reordered(bv).WeightedSpectraDKI=Reordered(bv).WeightedSpectraDKI;
    Reordered(bv).Direction(dir).nbkeeperDKI=0;
    end
    Reordered(bv).OutliersDKI.TotalKept=[Reordered(bv).OutliersDKI.TotalKept Reordered(bv).Direction(dir).nbkeeperDKI];
    end
 
    % computing new skewness and % of outliers
    Reordered(bv).WeightedSpectraDKI=Reordered(bv).WeightedSpectraDKI*nb_av/(nb_av-EjectedDirCounter);
    Reordered(bv).OutliersDKI.SkewnessSorted=skewness(Reordered(bv).OutliersDKI.SortedMeanSignal);
    Reordered(bv).OutliersDKI.Percent=1-sum(Reordered(bv).OutliersDKI.TotalKept)/size(Reordered(bv).Spectra,2);
end

% plotting mean upfield signal for all transients & threshold
figure
hold on 
scatter([1:size(ReorderedMean,2)],ReorderedMean, 10, 'k')
scatter([1:size(ReorderedMeanDownfield,2)],ReorderedMeanDownfield,10,'+','MarkerEdgeColor',[.7 .7 .7])
plot(Median_plot, 'Color',[0 0 0], 'Linewidth', 1.5)
plot(MAD_plot.B0, 'Color',[1 0.2 0.7], 'Linewidth', 1.2)
plot(MAD_plot.B0_plus,'--','Color',[1 0.2 0.7])
plot(MAD_plot.BV, 'Color',[0.5 0.2 0.4], 'Linewidth', 1.2)
plot(MAD_plot.BV_plus,'--','Color',[0.5 0.2 0.4])
plot(MAD_plot.DKI, 'Color',[0 1 1], 'Linewidth', 1.5)
plot(MAD_plot.DKI_plus,'--','Color',[0 1 1])


xlim([-50 620])    


fDKI=figure 
for bv=1:6
    subplot(2,6,bv)
    histogram(Reordered(bv).MeanSignal,'FaceColor',coloursHist(bv,:),'EdgeColor','none','BinWidth',FactorMAD*Reordered(1).MAD/20)
    subplot(2,6,bv+6)
    histogram(Reordered(bv).OutliersDKI.SortedMeanSignal,'FaceColor',coloursHist(bv,:),'EdgeColor','none','BinWidth',FactorMAD*Reordered(1).MAD/20)
end

name_fig=strcat(ProcessedFolder,filesep,'HistogramDKI.fig');
savefig(fDKI, name_fig)

%% 6D - Outlier detection based on the DKI prediction following 6A | median & MAD(bv=0)


MAD_plot.DKIBV=[];
MAD_plot.DKIBV_plus=[];
DKIBV_pred_plot=[];

for bv=1:6
Reordered(bv).WeightedSpectraDKIBV=[ ];
Reordered(bv).OutliersDKIBV.SortedMeanSignal=[ ];
Reordered(bv).OutliersDKIBV.TotalKept=[];
Reordered(bv).PredBV=[];
end


for bv=1
    MAD_plot.DKIBV=[MAD_plot.DKIBV (Reordered(bv).median_spec - FactorMAD*Reordered(1).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    MAD_plot.DKIBV_plus=[MAD_plot.DKIBV_plus (Reordered(bv).median_spec + FactorMAD*Reordered(1).MAD)*ones(1,Nb_Repetitions_per_b(bv)*nb_av*size(Acquisitions,2))];
    DKIBV_pred_plot=[DKI_pred_plot Reordered(bv).median_spec*ones(1,nb_av*Nb_Repetitions_per_b(bv)*size(Acquisitions,2))];
    for dir=1:nb_av
    Reordered(bv).Direction(dir).weightedspectraDKIBV=Reordered(bv).Direction(dir).weightedspectraB0;
    Reordered(bv).WeightedSpectraDKIBV=Reordered(bv).WeightedSpectraB0;
    
    Reordered(bv).Direction(dir).nbkeeperDKIBV=Reordered(bv).Direction(dir).nbkeeperB0;
    end
    Reordered(bv).OutliersDKIBV.SortedMeanSignal=Reordered(bv).OutliersB0.SortedMeanSignal;
    Reordered(bv).OutliersDKIBV.SkewnessSorted=Reordered(bv).OutliersB0.SkewnessSorted;
    Reordered(bv).OutliersDKIBV.Percent=Reordered(bv).OutliersB0.Percent;
  
end



for bv=2:6

    % sorting out at the direction & b-value levels (but threshold defined at the b-value level)

     EjectedDirCounter=0;
     
    for dir=1:nb_av
        
    MAD_plot.DKIBV=[MAD_plot.DKIBV (Reordered(bv).predBV(dir)*ones(1,Nb_Repetitions_per_b(bv)*2)+ Reordered(bv).median_spec - mean(Reordered(bv).predBV)    - FactorMAD*Reordered(1).MAD)];
    MAD_plot.DKIBV_plus=[MAD_plot.DKIBV_plus (Reordered(bv).predBV(dir)*ones(1,Nb_Repetitions_per_b(bv)*2)+ Reordered(bv).median_spec - mean(Reordered(bv).predBV)    + FactorMAD*Reordered(1).MAD)];
    DKIBV_pred_plot=[DKI_pred_plot Reordered(bv).predBV(dir)*ones(1,Nb_Repetitions_per_b(bv)*2)+ Reordered(bv).median_spec - mean(Reordered(bv).predBV)   ];
    
    Reordered(bv).Direction(dir).keepersDKIBV= find((Reordered(bv).Direction(dir).mean_signalBV -(Reordered(bv).predBV(dir)-mean(Reordered(bv).predBV)+Reordered(bv).median_spec)  + FactorMAD*Reordered(1).MAD)>0);
    

    
    if size(Reordered(bv).Direction(dir).keepersDKIBV,2)>0
    Reordered(bv).Direction(dir).weightedspectraDKIBV=Reordered(bv).Direction(dir).spectraBV(:,Reordered(bv).Direction(dir).keepersDKIBV)/size(Reordered(bv).Direction(dir).keepersDKIBV,2);
    Reordered(bv).WeightedSpectraDKIBV=[Reordered(bv).WeightedSpectraDKIBV Reordered(bv).Direction(dir).spectraBV(:,Reordered(bv).Direction(dir).keepersDKIBV)/size(Reordered(bv).Direction(dir).keepersDKIBV,2)];
    Reordered(bv).OutliersDKIBV.SortedMeanSignal=[Reordered(bv).OutliersDKIBV.SortedMeanSignal Reordered(bv).Direction(dir).mean_signalBV(Reordered(bv).Direction(dir).keepersDKIBV)];
    Reordered(bv).Direction(dir).nbkeeperDKIBV=size(Reordered(bv).Direction(dir).keepersDKIBV,2);
    else 
    EjectedDirCounter=EjectedDirCounter+1;
    Reordered(bv).WeightedSpectraDKIBV=Reordered(bv).WeightedSpectraDKIBV;
    Reordered(bv).Direction(dir).nbkeeperDKIBV=0;
    end
    Reordered(bv).OutliersDKIBV.TotalKept=[Reordered(bv).OutliersDKIBV.TotalKept Reordered(bv).Direction(dir).nbkeeperDKIBV];
    end
    Reordered(bv).WeightedSpectraDKIBV=Reordered(bv).WeightedSpectraDKIBV*nb_av/(nb_av-EjectedDirCounter);    % computing new skewness and % of outliers
    Reordered(bv).OutliersDKIBV.SkewnessSorted=skewness(Reordered(bv).OutliersDKIBV.SortedMeanSignal);
    Reordered(bv).OutliersDKIBV.Percent=1-sum(Reordered(bv).OutliersDKIBV.TotalKept)/size(Reordered(bv).Spectra,2);
end


% plotting mean upfield signal for all transients & threshold
scatterFig=figure
hold on 
scatter([1:size(ReorderedMean,2)],ReorderedMean, 10, 'k')
scatter([1:size(ReorderedMeanDownfield,2)],ReorderedMeanDownfield,10,'+','MarkerEdgeColor',[.7 .7 .7])
plot(Median_plot, 'Color',[0 0 0], 'Linewidth', 1.5)
plot(MAD_plot.B0, 'Color',[1 0.2 0.7], 'Linewidth', 1.2)
plot(MAD_plot.B0_plus,'--','Color',[1 0.2 0.7])
plot(MAD_plot.BV, 'Color',[0.5 0.2 0.4], 'Linewidth', 1.2)
plot(MAD_plot.BV_plus,'--','Color',[0.5 0.2 0.4])
plot(MAD_plot.DKI, 'Color',[0 .5 .5], 'Linewidth', 1.2)
plot(MAD_plot.DKI_plus,'--','Color',[0 .5 .5])
plot(MAD_plot.DKIBV, 'Color',[0 1 1], 'Linewidth', 1.5)
plot(MAD_plot.DKIBV_plus,'--','Color',[0 1 1])


xlim([-50 620])    

name_fig=strcat(ProcessedFolder,filesep,'Scatterfig.fig');
savefig(scatterFig, name_fig)


fDKIBV=figure 
for bv=1:6
    subplot(2,6,bv)
    histogram(Reordered(bv).MeanSignal,'FaceColor',coloursHist(bv,:),'EdgeColor','none','BinWidth',FactorMAD*Reordered(1).MAD/20)
    subplot(2,6,bv+6)
    histogram(Reordered(bv).OutliersDKIBV.SortedMeanSignal,'FaceColor',coloursHist(bv,:),'EdgeColor','none','BinWidth',FactorMAD*Reordered(1).MAD/20)
end

name_fig=strcat(ProcessedFolder,filesep,'HistogramDKIBV.fig');
savefig(fDKIBV, name_fig)


%% 7 - Summing spectra

% Summing spectra from 6A

for bv=1:6
    Reordered(bv).WeightedSpectraSingBV=[];
    for dir=1:nb_av
        if size(Reordered(bv).Direction(dir).weightedspectraBV,2) == 0
    Reordered(bv).WeightedSpectraSingBV=Reordered(bv).WeightedSpectraSingBV;
        else        

[~,fid_phasedtot,fid_rawtot,~]=phaser(Reordered(bv).Direction(dir).weightedspectraBV);
Reordered(bv).Direction(dir).meanspectraBV=fid_phasedtot;
Reordered(bv).WeightedSpectraSingBV=[Reordered(bv).WeightedSpectraSingBV fid_phasedtot'];
        end
    end

    [~,fid_phasedTOT,fid_rawtot,~]=phaser(Reordered(bv).WeightedSpectraSingBV);
    Reordered(bv).MeanSpectraBV=fid_phasedTOT;
    
end


figure
for bv=1:6
hold on 
plot(abs(fftshift(fft(Reordered(bv).MeanSpectraBV))));
end


% Summing spectra from 6B


for bv=1:6
 Reordered(bv).WeightedSpectraSingB0=[];   
    for dir=1:nb_av
        if size(Reordered(bv).Direction(dir).weightedspectraB0,2) == 0
    Reordered(bv).WeightedSpectraSingB0=Reordered(bv).WeightedSpectraSingB0;
        else

[~,fid_phasedtot,fid_rawtot,~]=phaser(Reordered(bv).Direction(dir).weightedspectraB0);
Reordered(bv).Direction(dir).meanspectraB0=fid_phasedtot;
Reordered(bv).WeightedSpectraSingB0=[Reordered(bv).WeightedSpectraSingB0 fid_phasedtot'];
        end
    end

    [~,fid_phasedTOT,fid_rawtot,~]=phaser(Reordered(bv).WeightedSpectraSingB0);
    Reordered(bv).MeanSpectraB0=fid_phasedTOT;
    
end


figure
for bv=1:6
hold on 
plot(abs(fftshift(fft(Reordered(bv).MeanSpectraB0))));
end


% Summing spectra from 6C

for bv=1:6
    
 Reordered(bv).WeightedSpectraSingDKI=[];   
    for dir=1:nb_av
        
        if size(Reordered(bv).Direction(dir).weightedspectraDKI,2) == 0
Reordered(bv).WeightedSpectraSingDKI=Reordered(bv).WeightedSpectraSingDKI;
        else
[~,fid_phasedtot,fid_rawtot,~]=phaser(Reordered(bv).Direction(dir).weightedspectraDKI);
Reordered(bv).Direction(dir).meanspectraDKI=fid_phasedtot;
Reordered(bv).WeightedSpectraSingDKI=[Reordered(bv).WeightedSpectraSingDKI fid_phasedtot'];
        end
    end

    [~,fid_phasedTOT,fid_rawtot,~]=phaser(Reordered(bv).WeightedSpectraSingDKI);
    Reordered(bv).MeanSpectraDKI=fid_phasedTOT;
    
end



figure
for bv=1:6
hold on 
plot(abs(fftshift(fft(Reordered(bv).MeanSpectraDKI))));
end

% Summing spectra from 6D

for bv=1:6
    
 Reordered(bv).WeightedSpectraSingDKIBV=[];   
    for dir=1:nb_av
        
        if size(Reordered(bv).Direction(dir).weightedspectraDKIBV,2) == 0
Reordered(bv).WeightedSpectraSingDKIBV=Reordered(bv).WeightedSpectraSingDKIBV;
        else
[~,fid_phasedtot,fid_rawtot,~]=phaser(Reordered(bv).Direction(dir).weightedspectraDKIBV);
Reordered(bv).Direction(dir).meanspectraDKIBV=fid_phasedtot;
Reordered(bv).WeightedSpectraSingDKIBV=[Reordered(bv).WeightedSpectraSingDKIBV fid_phasedtot'];
        end
    end

    [~,fid_phasedTOT,fid_rawtot,~]=phaser(Reordered(bv).WeightedSpectraSingDKI);
    Reordered(bv).MeanSpectraDKIBV=fid_phasedTOT;
    
end



figure
for bv=1:6
hold on 
plot(abs(fftshift(fft(Reordered(bv).MeanSpectraDKIBV))));
end



name_serie_acq=strcat(ProcessedFolder,filesep,'ReorderedData','.mat');
save(name_serie_acq,'Reordered', 'ReorderedMean', 'ReorderedMeanBV','ReorderedMeanDownfield','MAD_plot','Median_plot','DKI_pred_plot','DKIBV_pred_plot','-mat')


%% 8 - Save to LCModel prior to eddy current correction & water residual removal

figure
for bv=1:size(b_values,1)   

directory_LCModel_up=strcat('Fitting_LCModel');  

directory_LCModelBV = strcat(directory_LCModel_up,filesep,Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'BV',filesep);
mkdir(directory_LCModelBV)
save_fid(Reordered(bv).MeanSpectraBV',0,directory_LCModelBV);

directory_LCModelB0 = strcat(directory_LCModel_up,filesep,Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'B0',filesep);
mkdir(directory_LCModelB0)
save_fid(Reordered(bv).MeanSpectraB0',0,directory_LCModelB0);

directory_LCModelDKI = strcat(directory_LCModel_up,filesep,Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'DKI',filesep);
mkdir(directory_LCModelDKI)
save_fid(Reordered(bv).MeanSpectraDKI',0,directory_LCModelDKI);

directory_LCModelDKIBV = strcat(directory_LCModel_up,filesep,Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'DKIBV',filesep);
mkdir(directory_LCModelDKIBV)
save_fid(Reordered(bv).MeanSpectraDKIBV',0,directory_LCModelDKIBV);
   
end


%% 9 - Eddy current correction & save to LCModel post EC corr & water residual removal


try
    %%% Eddy Current correction & Water residual removal  
for bv=1:size(b_values,1)-2%k=5:6%
    [td_synth, ~, ~]=svdfid(Serie_water(bv,:)', 8, 6000, -150, 150, 60, 0.001, 2000);

    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraBV',td_synth);    
    [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    Reordered(bv).MeanSpectraBV_EC=td_diff';
    
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraB0',td_synth);    
    [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    Reordered(bv).MeanSpectraB0_EC=td_diff';
    
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraDKI',td_synth);    
    [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    Reordered(bv).MeanSpectraDKI_EC=td_diff';
    
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraDKIBV',td_synth);    
    [~, td_diff, ~]=svdfid(fid_cor, 8, 6000, -150, 150, 60, 0.001, 2000);
    Reordered(bv).MeanSpectraDKIBV_EC=td_diff';
   
end
    %%% Eddy Current correction only (water signal too low) 
for bv=5:6%
    [td_synth, ~, ~]=svdfid(Serie_water(bv,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
   
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraBV',td_synth);
    Reordered(bv).MeanSpectraBV_EC=fid_cor';
    
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraB0',td_synth);
    Reordered(bv).MeanSpectraB0_EC=fid_cor';
%     
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraDKI',td_synth);
    Reordered(bv).MeanSpectraDKI_EC=fid_cor';

    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraDKIBV',td_synth);
    Reordered(bv).MeanSpectraDKIBV_EC=fid_cor';
    
end

catch
    try    
    
    %%% Eddy Current correction only for all b (water signal too low)
for bv=1:6%
    [td_synth, ~, ~]=svdfid(Serie_water(bv,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraBV',td_synth);
    Reordered(bv).MeanSpectraBV_EC=fid_cor';
    
    [td_synth, ~, ~]=svdfid(Serie_water(bv,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraB0',td_synth);
    Reordered(bv).MeanSpectraB0_EC=fid_cor';
    
    [td_synth, ~, ~]=svdfid(Serie_water(bv,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraDKI',td_synth);
    Reordered(bv).MeanSpectraDKI_EC=fid_cor';

    [td_synth, ~, ~]=svdfid(Serie_water(bv,:)', 8, 6000, -150, 150, 60, 0.001, 2000);
    fid_cor=deconv_ECC(Reordered(bv).MeanSpectraDKIBV',td_synth);
    Reordered(bv).MeanSpectraDKIBV_EC=fid_cor';

   
end 
    catch
        %%% Only keeping non corrected spectra because EC corr failed. We
        %%% record in a 'EC_errors.txt' file which acquisitions could not
        %%% be corrected.
    Reordered(bv).MeanSpectraBV_EC=Reordered(bv).MeanSpectraBV;
    Reordered(bv).MeanSpectraB0_EC=Reordered(bv).MeanSpectraB0;
    Reordered(bv).MeanSpectraDKI_EC=Reordered(bv).MeanSpectraDKI;
    Reordered(bv).MeanSpectraDKIBV_EC=Reordered(bv).MeanSpectraDKIBV;
    line=strcat('Pup ', num2str(Pup(puppi)),', P',num2str(Age(agecounter)),', from litter ',num2str(Litter(litcounter)),' has not run Eddy current correction')
    filename=strcat(directory_LCModel_up,filesep,'EC_errors.txt');
   	fileID=fopen(filename,'a')
    fprintf(fileID,'%s\n',line)
    fclose(fileID)
    end

end


%%% Save to LCModel the corrected acquisitions
for bv=1:size(b_values,1)   
  
directory_LCModelBV_EC = strcat(directory_LCModel_up,filesep,'EC_',Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'BV',filesep);
mkdir(directory_LCModelBV_EC)
save_fid(Reordered(bv).MeanSpectraBV_EC',0,directory_LCModelBV_EC);

directory_LCModelB0_EC = strcat(directory_LCModel_up,filesep,'EC_',Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'B0',filesep);
mkdir(directory_LCModelB0_EC)
save_fid(Reordered(bv).MeanSpectraB0_EC',0,directory_LCModelB0_EC);

directory_LCModelDKI_EC = strcat(directory_LCModel_up,filesep,'EC_',Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'DKI',filesep);
mkdir(directory_LCModelDKI_EC)
save_fid(Reordered(bv).MeanSpectraDKI_EC',0,directory_LCModelDKI_EC);

directory_LCModelDKIBV_EC = strcat(directory_LCModel_up,filesep,'EC_',Name_Acquisition,filesep,'TM',num2str(protocol(Acquisitions(1)).TM),filesep,'b',num2str(b_values(bv)),filesep,'DKIBV',filesep);
mkdir(directory_LCModelDKIBV_EC)
save_fid(Reordered(bv).MeanSpectraDKIBV_EC',0,directory_LCModelDKIBV_EC);
   

end



name_serie_acq=strcat(ProcessedFolder,filesep,'ReorderedData','.mat');
save(name_serie_acq,'Reordered', 'ReorderedMean', 'ReorderedMeanBV','ReorderedMeanDownfield','MAD_plot','Median_plot','DKI_pred_plot','DKIBV_pred_plot','-mat')
