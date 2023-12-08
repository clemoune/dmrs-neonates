%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading Bruker Protocols and Running Processing for High b-values OR Long
% TM acquired in rat neonates.
% Based on the read_parameters function written by Matteo Caffini (github)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DataType is 'highb' OR 'longtm'
DataType='highb' 
addpath('support_functions/')
data_nb_dir='/vols/Data/preclinical/Clemence/Purkinje_SeedGrant/';
Litter=[1 2 3];
Age=[5 10 15 20 30];
Pup=[1:7];

for litcounter=1:size(Litter,2)
    
    for agecounter=1:size(Age,2)

        for puppi=1:size(Pup,2)
    
    close all
    clearvars -except data_nb_dir Pup puppi Litter litcounter Age agecounter
    clc

    tic    
    
Data_Dir=strcat(data_nb_dir,filesep,'Litter',num2str(Litter(litcounter)),'_Rat',filesep,'Neonates_P',num2str(Age(agecounter)),filesep,'P',num2str(Age(agecounter)),'_Pup',num2str(Pup(puppi)));

            if exist(Data_Dir, 'dir') == 7

Number_Acq_max=50;

strcat('Pup ', num2str(Pup(puppi)),', P',num2str(Age(agecounter)),', from litter ',num2str(Litter(litcounter)),' is starting...')

for i=1:Number_Acq_max
    
    puper_name=strcat(Data_Dir,filesep,num2str(i));

    if exist(puper_name)==7
        parameters_method=read_parameters(strcat(puper_name,filesep,'method'));
        protocol(i).type=parameters_method.Method;
        protocol(i).TE=parameters_method.PVM_EchoTime;
        protocol(i).TR=parameters_method.PVM_RepetitionTime;
        protocol(i).Nb_Averages=parameters_method.PVM_NAverages;
        protocol(i).Nb_Repetitions=parameters_method.PVM_NRepetitions;
        protocol(i).Nb_Repetitions=parameters_method.PVM_NRepetitions;
        if matches(parameters_method.Method,'User:cl_STELASER_PA360_b')==1
            protocol(i).b_value=parameters_method.B_values_list;
            protocol(i).rep_per_b=parameters_method.Repetitions_per_b_list;
            protocol(i).TE_STE=parameters_method.TE_STE;
            protocol(i).TM=parameters_method.MixingTime;
            protocol(i).Macro=parameters_method.DIR_Module;
            protocol(i).WaterSup=parameters_method.PVM_WsMode;

        end
        
        if matches(parameters_method.Method,'User:cl_STELASER_PA360')==1
            protocol(i).b_value=parameters_method.B_value;
            %protocol(i).rep_per_b=parameters_method.Repetitions_per_b_list;
            protocol(i).TE_STE=parameters_method.TE_STE;
            protocol(i).TM=parameters_method.MixingTime;
            protocol(i).Macro=parameters_method.DIR_Module;
            protocol(i).WaterSup=parameters_method.PVM_WsMode;

        end        
    else
    end
     

end

%%% Running the processing (calling scripts)

Name_Acquisition = strcat('Lit',num2str(Litter(litcounter)),'P',num2str(Age(agecounter)),'Pup',num2str(Pup(puppi)));

if DataType=='highb'
BSB_HighB_Pups
else
BSB_PhaseEstimated
BSB_LongTM_Pups
end


time_elapsed=toc

strcat('Pup ', num2str(Pup(puppi)),', P',num2str(Age(agecounter)),', from litter ',num2str(Litter(litcounter)),' is processed!')
            end 
        end
    end 
end

