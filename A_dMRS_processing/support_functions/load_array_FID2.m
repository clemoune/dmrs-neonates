function [FID]=load_array_FID2(File_Name,Nfid);

% nb_pts_to_remove=68;
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

%cd /neurospin/iseult/julien_V

%cd('D:\Documents and Settings\jv201216\Bureau\CEA');

%[FileName,PathName] = uigetfile('*.*',text);
%signal_input = fullfile( PathName, FileName ) ;
openMRI = fopen(File_Name,'r');

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

np=length(FID)/Nfid;

FID=reshape(FID,[np Nfid]);

FID=[FID(nb_pts_to_remove+1:end,:); zeros(nb_pts_to_remove,Nfid)];

cd(currentdir);