function [SS_pred, sd_SS_pred]=plot_fit_ClemenceCerebellum(metabolite_meas, exp_T,fit_result,dictionary, plot_flag, plot_color)

%--------------------------------------------------------------------------
%
% Plot the Random Forest estimated best fit curve and output the curve in
% SS_pred.
%
% metabolite_meas: the experimental data corresponding to the desired
% metabolite. 
%
% fit_result: the vectore of 5 parameters estimated by Random Forest
% regression. For example, [Dintra N L SDN SDL]
%
% plot_flag: 1 to let the function plot the estimated curve fitted to the
% experimental data
%
% [SS_pred, sd_SS_pred]: predicted ADC as a function of diffusion time, and
% corresponding uncertainty, respectively.
%
%--------------------------------------------------------------------------

% load Dictionary.mat
% 
% Y_full{1}=Y_dictionary;
% Dintra_full{1}=Dintra;
% N_full{1}=N;
% L_full{1}=L;
% SDN_full{1}=SDN;
% SDL_full{1}=SDL;
% 
% load Dictionary_2_partial.mat
% 
% Y_full{2}=Y_dictionary;
% Dintra_full{2}=Dintra;
% N_full{2}=N;
% L_full{2}=L;
% SDN_full{2}=SDN;
% SDL_full{2}=SDL;
% 
% load Dictionary_3.mat
% 
% Y_full{3}=Y_dictionary;
% Dintra_full{3}=Dintra;
% N_full{3}=N;
% L_full{3}=L;
% SDN_full{3}=SDN;
% SDL_full{3}=SDL;
% 
% load Dictionary_Part_II.mat
% 
% Y_full{4}=Y_dictionary;
% Dintra_full{4}=Dintra;
% N_full{4}=N;
% L_full{4}=L;
% SDN_full{4}=SDN;
% SDL_full{4}=SDL;
% 
% load Dictionary_Part_JV.mat
% 
% Y_full{5}=Y_dictionary;
% Dintra_full{5}=Dintra;
% N_full{5}=N;
% L_full{5}=L;
% SDN_full{5}=SDN;
% SDL_full{5}=SDL;

Y_full = dictionary.Y_full;
Dintra_full = dictionary.Dintra_full;
N_full = dictionary.N_full;
L_full = dictionary.L_full;
SDN_full = dictionary.SDN_full;
SDL_full = dictionary.SDL_full;

Signals_FV = [];

Params = [];

kk=1;

for dd = 1:5
    for i=1:length(Dintra_full{dd})
        for j=1:length(N_full{dd})
            for k=1:length(L_full{dd})
                for l=1:length(SDN_full{dd})
                    for m=1:length(SDL_full{dd})
                        
                        Y_tmp(1:2)=Y_full{dd}([2 5],i,j,k,l,m);
                        Y_tmp(3)=0.7.*Y_full{dd}(6,i,j,k,l,m) + 0.3.*Y_full{dd}(7,i,j,k,l,m);
                        Y_tmp(4)=Y_full{dd}(7,i,j,k,l,m);
                        if sum(Y_tmp)>0
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

Y_exp = metabolite_meas;

rr = fit_result;

YY=Params-rr;
Obj_fun = sqrt(sum(YY.^2,2));

[~, I] = sort(Obj_fun);

SS_reord = Signals_FV(I,:);

SS_pred = nanmean(SS_reord(1:4,:),1);
sd_SS_pred = zeros(1,size(Signals_FV,2));

for i=1:size(Signals_FV,2)
    
    sd_SS_pred(i) = nanstd(SS_reord(1:4,i))/2;
    
end

if plot_flag==1
    hold on
    plot(exp_T, Y_exp, 'ko', 'markersize', 12, 'markerfacecolor','k')
    plot(exp_T, SS_pred, plot_color, 'linewidth',2)
    %plot(exp_T, SS_pred+sd_SS_pred, 'b--', 'linewidth',3)
    %plot(exp_T, SS_pred-sd_SS_pred, 'b--', 'linewidth',3)
    axis([0 1200 0 0.3])
end

end
