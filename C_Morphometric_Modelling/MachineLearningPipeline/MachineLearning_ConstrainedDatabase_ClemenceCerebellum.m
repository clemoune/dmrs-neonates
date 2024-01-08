function []=MachineLearning_ConstrainedDatabase_ClemenceCerebellum(Dintra_bounds, N_bounds, L_bounds,SDN_bounds,SDL_bounds)

load Dictionary.mat

Y_full{1}=Y_dictionary;
Dintra_full{1}=Dintra;
N_full{1}=N;
L_full{1}=L;
SDN_full{1}=SDN;
SDL_full{1}=SDL;

load Dictionary_2_partial.mat

Y_full{2}=Y_dictionary;
Dintra_full{2}=Dintra;
N_full{2}=N;
L_full{2}=L;
SDN_full{2}=SDN;
SDL_full{2}=SDL;

load Dictionary_3.mat

Y_full{3}=Y_dictionary;
Dintra_full{3}=Dintra;
N_full{3}=N;
L_full{3}=L;
SDN_full{3}=SDN;
SDL_full{3}=SDL;

load Dictionary_Part_II.mat

Y_full{4}=Y_dictionary;
Dintra_full{4}=Dintra;
N_full{4}=N;
L_full{4}=L;
SDN_full{4}=SDN;
SDL_full{4}=SDL;

load Dictionary_Part_JV.mat

Y_full{5}=Y_dictionary;
Dintra_full{5}=Dintra;
N_full{5}=N;
L_full{5}=L;
SDN_full{5}=SDN;
SDL_full{5}=SDL;


close all

%  Y_dictionary: [8    10   7  10    4    5]
%                [td Dintra N   L   SDN  SDL]

% num_curves = size(Y_dictionary,2)*size(Y_dictionary,3)*size(Y_dictionary,4)*size(Y_dictionary,5)*size(Y_dictionary,6);
% num_test=2500;
% num_train=num_curves-num_test;

Signals_FV = [];
% TrainSignals_FV = zeros(num_train, 6);
% TestSignals_FV = zeros(num_test, 6);

Params = [];
% TrainParams = zeros(num_train, 5);
% TestParams = zeros(num_test, 5);

kk=1;

for dd = 1:5
    for i=1:length(Dintra_full{dd})
        for j=1:length(N_full{dd})
            for k=1:length(L_full{dd})
                for l=1:length(SDN_full{dd})
                    for m=1:length(SDL_full{dd})
                        
                        Y_tmp(1:2)=Y_full{dd}([2 5],i,j,k,l,m);
                        Y_tmp(3)= 0.7.*Y_full{dd}(6,i,j,k,l,m) + 0.3.*Y_full{dd}(7,i,j,k,l,m);
                        Y_tmp(4)=Y_full{dd}(7,i,j,k,l,m);
                        
                        Dintra_condition = Dintra_full{dd}(i)>=Dintra_bounds(1) & Dintra_full{dd}(i)<=Dintra_bounds(2);
                        N_condition = N_full{dd}(j)>=N_bounds(1) & N_full{dd}(j)<=N_bounds(2);
                        L_condition = L_full{dd}(k)>=L_bounds(1) & L_full{dd}(k)<=L_bounds(2);
                        SDN_condition = SDN_full{dd}(l)>=SDN_bounds(1) & SDN_full{dd}(l)<=SDN_bounds(2);
                        SDL_condition = SDL_full{dd}(m)>=SDL_bounds(1) & SDL_full{dd}(m)<=SDL_bounds(2);
                        
                        total_condition = Dintra_condition & N_condition & L_condition & SDN_condition & SDL_condition;
                        
                        if ~isempty(Y_tmp) && sum(Y_tmp)>0
                            if Dintra_condition && N_condition && L_condition && SDN_condition && SDL_condition
                                Signals_FV(kk,:)=Y_tmp;
                                Params(kk,:)=[Dintra_full{dd}(i) N_full{dd}(j) L_full{dd}(k) SDN_full{dd}(l) SDL_full{dd}(m)];
                                kk=kk+1;
                            end
                        end
                    end
                end
            end
        end
    end
end
num_curves = size(Signals_FV,1);
num_test=2500;
num_train=num_curves-num_test;

TrainSignals_FV = zeros(num_train, length(Y_tmp));
% TestSignals_FV = zeros(num_test, 6);

TrainParams = zeros(num_train, 5);
% TestParams = zeros(num_test, 5);

II=randi(num_curves,[num_test 1]);

TestSignals_FV=Signals_FV(II,:);
TestParams=Params(II,:);

JJ=(1:num_curves)';

kk=1;

for i=1:num_curves
    if sum(JJ(i)==II)==0
        
        TrainSignals_FV(kk,:)=Signals_FV(i,:);
        TrainParams(kk,:)=Params(i,:);
        kk=kk+1;
    end
end

save TrainSignals_FV_extended_ClemenceCerebellum_constrained TrainSignals_FV
save TrainParams_extended_ClemenceCerebellum_constrained TrainParams

save TestSignals_FV_extended_ClemenceCerebellum_constrained TestSignals_FV
save TestParams_extended_ClemenceCerebellum_constrained TestParams

RF_features_reduced=Signals_FV;
param_to_fit = Params;

save RF_features_extended_ClemenceCerebellum_constrained RF_features_reduced
save param_to_fit_extended_ClemenceCerebellum_constrained param_to_fit

end
