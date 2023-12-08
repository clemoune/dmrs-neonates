function [FID1 FID2 FID3 FID4]=coil_decomposition(FileName,nb_points,nb_repetitions,nb_coils);

% FileName = complete path of your data 'C\...fullpath\file_acq\rawdata.job0'
% nb_points = number of points in your FID
% nb_repetitions = number of repetitions
% nb_coils = number of coils

%nb_pts_to_remove=68;
nb_pts_to_remove=76;

%******************************************
% loadMRI.m
% 
% Reads BRUKER FID file and converts it to
% MATLAB PC format.
%******************************************

%*******************************************************
% Read BRUKER file (already byte-swapped with pcproc)
%*******************************************************
%%% Open file dialog box %%

currentdir=pwd;

openMRI = fopen(FileName,'r');

if (openMRI < 0)
   figure('Position',[400 400 400 100]);

   axes('Visible','off');

   text('Position',[0.0, 0.60],'String','File does not exist or can not be opened !!!');
else
    indata1 = fread(openMRI, inf, 'long');
    fclose(openMRI);

    %****************************************************
    % Make (real,imag) pairs and store in complex matrix
    %****************************************************
    size_FID = size(indata1,1)/2;
    temp1 = reshape(indata1,2,size_FID);
    temp2 = temp1(1,:) + 1i*temp1(2,:);
    FID = reshape(temp2, size_FID, 1);
    size(FID);

    clear temp1 temp2;

end

FID=[FID(nb_pts_to_remove+1:end); zeros(nb_pts_to_remove,1)];


Divided_FID=zeros;
total_size=size(FID);
nb_FID=total_size(1)/nb_points;

Divided_FID_coil_aver=zeros;


FID1=zeros;
FID2=zeros;
FID3=zeros;
FID4=zeros;


for i=1:1:nb_FID
Divided_FID(1:nb_points,i)=FID((i-1)*nb_points+1:i*nb_points,1);
end 

for j=1:1:nb_repetitions;
FID1(1:nb_points,j)=Divided_FID(:,nb_coils*(j-1)+1);
FID2(1:nb_points,j)=Divided_FID(:,nb_coils*(j-1)+2);
FID3(1:nb_points,j)=Divided_FID(:,nb_coils*(j-1)+3);
FID4(1:nb_points,j)=Divided_FID(:,nb_coils*(j-1)+4);
end

cd(currentdir);