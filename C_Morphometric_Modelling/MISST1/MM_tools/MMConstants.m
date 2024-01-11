function protocol=MMConstants(model,protocol)
%
% camino.m--------------------------------------------------------------
% Calculates constants needed for the Matrix Method
% 
% protocol=MMConstants(model,protocol)
% 
% Description: Adds to the diffusion protocol the fields necessary for 
% evaluating the diffusion signal using the Matrix Method formalism. 
% Such fields include the matrices S,A and R which depend on the restriction
% specified in model.
%
% Parameters:
% protocol - structure describing the diffusion protocol  
% model - structure which describes the diffusion substrate 
%       model.name - model name which includes the different compartments
%       model.params - parameter values for the given model
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Ivana Drobnjak (i.drobnjak@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
% 

% Setting qunit, N, tau, gstep
GAMMA = 2.675987E8;

if ~isfield(protocol,'dimA')
  protocol.dimA=20;
end
if ~isfield(protocol,'diff')
  protocol.diff=1; % calculates the derivatives of the matrices
end
if ~isfield(protocol,'tau')
    if strcmp(protocol.pulseseq,'GEN') % needed when optimising each point on the waveform;
        protocol.tau=protocol.riset*(2+0.1);
    else
        protocol.tau=1E-4;
    end
end
if ~isfield(protocol,'mirror')
  protocol.mirror=0; % when 1 the gradient waveform is mirrored about the 180 degree pulse
end
if ~isfield(protocol,'complex')
    protocol.complex='abs'; % 
end
if isfield(protocol,'smalldel') && isfield(protocol,'delta') 
protocol.K = floor((protocol.smalldel+1E-10)./protocol.tau)+1;
protocol.dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
end
protocol.gstep=0.04/100;
protocol.qunit=GAMMA*protocol.tau*protocol.gstep;
dimA = protocol.dimA;
qunit=protocol.qunit/(2*pi);
protocol.pert=1E-5;
protocol.approx = 'MM';

strings = GetParameterStrings(model.name);

if sum(strcmp('rad',strings))  % cylinder restriction
    if sum(strcmp('shape',strings)) == 0 && sum(strcmp('var',strings)) == 0  % one radius
        radind = GetParameterIndex(model.name,'rad');
        [urad_vec urad_ind] = unique(model.params(:,radind), 'stable');
        modelcode = 1; % cylinder
        for j = 1:length(urad_vec)
            D = model.params(urad_ind(j),GetParameterIndex(model.name,'di')); % diffusivity
            a = model.params(urad_ind(j),GetParameterIndex(model.name,'rad')); % rad
            roo=MMroots(dimA,a,modelcode); % cylinder
            [protocol.S_cyl{j} protocol.A_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,a,dimA,roo,modelcode);
            protocol.R_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
            if protocol.diff==1
                % Radius apert
                apert=a*(1+protocol.pert);
                roo=MMroots(dimA,apert,modelcode);
                [protocol.Sperta_cyl{j} protocol.Aperta_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,apert,dimA,roo,modelcode);
                protocol.Rperta_cyl{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D); 
                % Dpert
                Dpert=D*(1+protocol.pert);
                roo=MMroots(dimA,a,modelcode);
                protocol.RpertD_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
            end
        end
    else
        radind = GetParameterIndex(model.name,'rad');
        shapeind = GetParameterIndex(model.name,'shape');
        varind = GetParameterIndex(model.name,'var');      
            if shapeind > 0
                shape = model.params(shapeind);
                b = model.params(radind)./shape;
                if ~isfield(protocol,'fixedR')
                     if ~isfield(protocol,'GammaRdist') || strcmp(protocol.GammaRdist,'pdf') % weights come from the pdf of gamma function
                        intR_steps = 10;
                        [R_vec, weight] = GammaRadList(shape, b, intR_steps);
                     elseif strcmp(protocol.GammaRdist,'vol') % diffusion signal from different radii is weighted by volume
                        intR_steps = 10;
                        if sum(strcmp('ecc',strings))
                            [R_vec, weight0] = GammaRadList(shape, b, intR_steps);
                            weight = weight0.*(R_vec*1E6).^3;
                            weight = weight./sum(weight);
                        else
                            [R_vec, weight] = GammaRadByVolList(shape, b, intR_steps);
                        end
                     else error Unknown GammaRdist option
                     end
                else
                     R_vec = protocol.GammaR;
                     weight = protocol.GammaWeights;   
                end
                 protocol.Rweight = weight;            
                if protocol.diff==1
                    bperta =  b*(1+protocol.pert);
                    bpertsh = b./(1+protocol.pert);
                    shapepert = shape*(1+protocol.pert);
                    if ~isfield(protocol,'GammaRdist') || strcmp(protocol.GammaRdist,'pdf') % weights come from the pdf of gamma function
                        intR_steps = 10;                        
                        [R_vecperta, weight_perta] = GammaRadList(shape, bperta, intR_steps);
                        [R_vecpertsh, weight_pertsh] = GammaRadList(shapepert, bpertsh, intR_steps);
                     elseif strcmp(protocol.GammaRdist,'vol') % diffusion signal from different radii is weighted by volume
                        intR_steps = 10;
                        [R_vecperta, weight_perta] = GammaRadList(shape, bperta, intR_steps);
                        [R_vecpertsh, weight_pertsh] = GammaRadList(shapepert, bpertsh, intR_steps);
                     else error Unknown GammaRdist option
                    end
                     protocol.Rweight_perta = weight_perta; 
                     protocol.Rweight_pertsh = weight_pertsh; 
                    
                end
            end
            if varind >0 
                if ~isfield(protocol,'fixedR')
                     if ~isfield(protocol,'NormalRdist') || strcmp(protocol.NormalRdist,'pdf') % weights come from the pdf of gamma function
                        intR_steps = 10;
                        [R_vec, weight] = NormalRadList(model.params(radind), model.params(varind), intR_steps);
                     elseif strcmp(protocol.NormalRdist,'vol') % diffusion signal from different radii is weighted by volume
                        intR_steps = 10;
                        [R_vec, weight] = NormalRadByVolList(model.params(radind),model.params(varind), intR_steps);
                     else error Unknown NormalRdist option

                    end
                else
                     R_vec = protocol.NormalR;
                     weight = protocol.NormalWeights;                     
                end
                protocol.Rweight = weight;               
                 if protocol.diff==1
                    varpert = model.params(varind)*(1+protocol.pert);
                    rpert = model.params(radind).*(1+protocol.pert);                    
                    if ~isfield(protocol,'GammaRdist') || strcmp(protocol.GammaRdist,'pdf') % weights come from the pdf of gamma function
                        intR_steps = 10;                        
                        [R_vecperta, weight_perta] =  NormalRadList(rpert, model.params(varind), intR_steps);
                        [R_vecpertvar, weight_pertvar] =  NormalRadList(model.params(radind), varpert, intR_steps);
                     elseif strcmp(protocol.GammaRdist,'vol') % diffusion signal from different radii is weighted by volume
                        intR_steps = 10;
                        [R_vecperta, weight_perta] = NormalRadByVolList(rpert, model.params(varind), intR_steps);
                        [R_vecpertvar, weight_pertvar] = NormalRadByVolList(model.params(radind), varpert, intR_steps);
                     else error Unknown GammaRdist option
                    end
                     protocol.Rweight_perta = weight_perta; 
                     protocol.Rweight_pertvar = weight_pertvar; 
                    
                end

            end
            
            
        modelcode = 1; % cylinder
        D = model.params(GetParameterIndex(model.name,'di')); % diffusivity
        for j = 1:length(R_vec)          
            a = R_vec(j); % rad
            roo=MMroots(dimA,a,modelcode); % cylinder
            [protocol.S_cyl{j} protocol.A_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,a,dimA,roo,modelcode);
            protocol.R_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
            if protocol.diff==1
                % mean Radius apert
                aperta=R_vecperta(j);
                roo=MMroots(dimA,aperta,modelcode);
                [protocol.Sperta_cyl{j} protocol.Aperta_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,aperta,dimA,roo,modelcode);
                protocol.Rperta_cyl{j}=MatrixR(protocol.tau,aperta,dimA,roo,modelcode,D); 
                % Dpert
                Dpert=D*(1+protocol.pert);
                roo=MMroots(dimA,a,modelcode);
                protocol.RpertD_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
                if shapeind > 0 % shape
                    apertsh=R_vecpertsh(j);
                    roo=MMroots(dimA,apertsh,modelcode);
                    [protocol.Spertsh_cyl{j} protocol.Apertsh_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,apertsh,dimA,roo,modelcode);
                    protocol.Rpertsh_cyl{j}=MatrixR(protocol.tau,apertsh,dimA,roo,modelcode,D); 
                    
                end
                 if varind > 0 % variance
                    apertvar=R_vecpertvar(j);
                    roo=MMroots(dimA,apertvar,modelcode);
                    [protocol.Spertvar_cyl{j} protocol.Apertvar_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,apertvar,dimA,roo,modelcode);
                    protocol.Rpertvar_cyl{j}=MatrixR(protocol.tau,apertvar,dimA,roo,modelcode,D); 
                    
                end
            end
        end
    end
end
if sum(strcmp('ecc',strings)) && sum(strcmp('lx',strings)) == 0  % planar restriction fot finite cylinders
    modelcode = 2; % planar
    if sum(strcmp('shape',strings)) == 0 && sum(strcmp('var',strings)) == 0  % one radius
        for j = 1:length(urad_vec)
            D = model.params(urad_ind(j),GetParameterIndex(model.name,'di')); % diffusivity
            a = model.params(urad_ind(j),GetParameterIndex(model.name,'rad')); % rad
            ecc = model.params(urad_ind(j),GetParameterIndex(model.name,'ecc')); % half the length; in the code the parallel planes are situated between -l and l
            l = a*ecc;
            roo=MMroots(dimA,l,modelcode); % planar
            [protocol.S_plane{j} protocol.A_plane{j} ]=MatrixSA(qunit,l,dimA,roo,modelcode);
            protocol.R_plane{j}=MatrixR(protocol.tau,l,dimA,roo,modelcode,D);
            if protocol.diff==1
                % Dpert
                Dpert = D*(1+protocol.pert);
                roo=MMroots(dimA,l,modelcode);
                protocol.RpertD_plane{j}=MatrixR(protocol.tau,l,dimA,roo,modelcode,Dpert); 
                % a pert
                lpert=l*(1+protocol.pert); % ecc * a*(1+pert)
                roo=MMroots(dimA,lpert,modelcode);
                [protocol.Sperta_plane{j} protocol.Aperta_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,modelcode);
                protocol.Rperta_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,modelcode,D); 
                % Ecc pert
                lpert=l*(1+protocol.pert); % ecc*(1+pert) * a
                roo=MMroots(dimA,lpert,modelcode);
                [protocol.SpertEcc_plane{j} protocol.ApertEcc_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,modelcode);
                protocol.RpertEcc_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,modelcode,D); 
            end
        end
    else
      
        % R_vec and the weights have been previously calculated; no need to repeat
         for j = 1:length(R_vec)
            D = model.params(GetParameterIndex(model.name,'di')); % diffusivity
            a = R_vec(j); % rad
            ecc = model.params(GetParameterIndex(model.name,'ecc')); % half the length; in the code the parallel planes are situated between -l and l
            l = a*ecc;
            roo=MMroots(dimA,l,modelcode); % planar
            [protocol.S_plane{j} protocol.A_plane{j} ]=MatrixSA(qunit,l,dimA,roo,modelcode);
            protocol.R_plane{j}=MatrixR(protocol.tau,l,dimA,roo,modelcode,D);
            if protocol.diff==1
                % Dpert
                Dpert = D*(1+protocol.pert);
                roo=MMroots(dimA,l,modelcode);
                protocol.RpertD_plane{j}=MatrixR(protocol.tau,l,dimA,roo,modelcode,Dpert); 
                % a pert
                lpert=ecc*R_vecperta(j); % ecc * a*(1+pert)
                roo=MMroots(dimA,lpert,modelcode);
                [protocol.Sperta_plane{j} protocol.Aperta_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,modelcode);
                protocol.Rperta_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,modelcode,D); 
                % Ecc pert
                lpert=l*(1+protocol.pert); % ecc*(1+pert) * a
                roo=MMroots(dimA,lpert,modelcode);
                [protocol.SpertEcc_plane{j} protocol.ApertEcc_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,modelcode);
                protocol.RpertEcc_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,modelcode,D); 
                 if shapeind > 0 % shape
                    lpert=R_vecpertsh(j)*ecc;
                     roo=MMroots(dimA,lpert,modelcode);
                    [protocol.Spertsh_plane{j} protocol.Apertsh_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,modelcode);
                    protocol.Rpertsh_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,modelcode,D); 

                  end
                 if varind > 0 % variance
                   lpert=R_vecpertvar(j)*ecc;
                     roo=MMroots(dimA,lpert,modelcode);
                    [protocol.Spertvar_plane{j} protocol.Apertvar_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,modelcode);
                    protocol.Rpertvar_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,modelcode,D); 
                    
                end
            
            end     
        end
        
        
    end
end
if sum(strcmp('rads',strings))  % spherical restriction
    radind = GetParameterIndex(model.name,'rads');
    [urad_vec urad_ind] = unique(model.params(:,radind), 'stable');
    modelcode = 0; % cylinder
    for j = 1:length(urad_vec)
        D = model.params(urad_ind(j),GetParameterIndex(model.name,'di')); % diffusivity
        a = model.params(urad_ind(j),GetParameterIndex(model.name,'rads')); % rad
        roo=MMroots(dimA,a,modelcode); % cylinder
        [Stemp Atemp protocol.kDn_sph]=MatrixSA(qunit,a,dimA,roo,modelcode);
        protocol.A_sph{j,1} = Atemp{1}; % A0
        protocol.A_sph{j,2} = Atemp{2}; % A90
        protocol.S_sph{j,1} = Stemp{1}; % S0
        protocol.S_sph{j,2} = Stemp{2}; % S90
        protocol.R_sph{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
        if protocol.diff==1
            % Radius apert
            apert=a*(1+protocol.pert);
            roo=MMroots(dimA,apert,modelcode);
            [Stemp Atemp ]=MatrixSA(qunit,apert,dimA,roo,modelcode);
            protocol.Aperta_sph{j,1} = Atemp{1}; % A0
            protocol.Aperta_sph{j,2} = Atemp{2}; % A90
            protocol.Sperta_sph{j,1} = Stemp{1}; % S0
            protocol.Sperta_sph{j,2} = Stemp{2}; % S90
            protocol.Rperta_sph{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D); 
            % Dpert
            Dpert=D*(1+protocol.pert);
            roo=MMroots(dimA,a,modelcode);
            protocol.RpertD_sph{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
        end
    end
end  
if (sum(strcmp('lx',strings))+sum(strcmp('ly',strings))+sum(strcmp('lz',strings)))==3  % cuboids
    lxind = GetParameterIndex(model.name,'lx');
    lyind = GetParameterIndex(model.name,'ly');
    lzind = GetParameterIndex(model.name,'lz');
    [ul_vec ul_ind] = unique([model.params(:,lxind) model.params(:,lyind) model.params(:,lzind)],'rows','stable');
    modelcode = 2; % planar
    for j = 1:length(ul_ind)        
        D = model.params(ul_ind(j),GetParameterIndex(model.name,'di')); % diffusivity
        a = model.params(ul_ind(j),GetParameterIndex(model.name,'lx'))/2;  % lx/2
        b = model.params(ul_ind(j),GetParameterIndex(model.name,'ly'))/2;  % ly/2
        c = model.params(ul_ind(j),GetParameterIndex(model.name,'lz'))/2;  % lz/2
        roo=MMroots(dimA,a,modelcode); % x dir
        [protocol.S_x{j} protocol.A_x{j} ]=MatrixSA(qunit,a,dimA,roo,modelcode);
        roo = MMroots(dimA,b,modelcode); % y dir
        [protocol.S_y{j} protocol.A_y{j} ]=MatrixSA(qunit,b,dimA,roo,modelcode);
        roo = MMroots(dimA,c,modelcode); % z dir
        [protocol.S_z{j} protocol.A_z{j} ]=MatrixSA(qunit,c,dimA,roo,modelcode);   
        % from roo MatrixR takes only the eigenvalues which are the same
        % no mater what length
        protocol.R_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
        protocol.R_y{j}=MatrixR(protocol.tau,b,dimA,roo,modelcode,D);
        protocol.R_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,D);
        if protocol.diff==1
            % lx pert
            apert=a*(1+protocol.pert);
            roo=MMroots(dimA,apert,modelcode);
            [protocol.Sperta_x{j} protocol.Aperta_x{j} ]=MatrixSA(qunit,apert,dimA,roo,modelcode);
            % ly pert
            bpert=b*(1+protocol.pert);
            roo=MMroots(dimA,bpert,modelcode);
            [protocol.Spertb_y{j} protocol.Apertb_y{j} ]=MatrixSA(qunit,bpert,dimA,roo,modelcode);
            %lz pert
            cpert=c*(1+protocol.pert);
            roo=MMroots(dimA,cpert,modelcode);
            [protocol.Spertc_z{j} protocol.Apertc_z{j} ]=MatrixSA(qunit,cpert,dimA,roo,modelcode);
      
            protocol.Rperta_x{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D); 
            protocol.Rpertb_y{j}=MatrixR(protocol.tau,bpert,dimA,roo,modelcode,D); 
            protocol.Rpertc_z{j}=MatrixR(protocol.tau,cpert,dimA,roo,modelcode,D); 
            % Dpert
            Dpert=D*(1+protocol.pert);
            protocol.RpertD_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
            protocol.RpertD_y{j}=MatrixR(protocol.tau,b,dimA,roo,modelcode,Dpert); 
            protocol.RpertD_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,Dpert); 
        end
    end
end
if (sum(strcmp('lx',strings))+sum(strcmp('ecc',strings))) == 2  % cuboids with lx = ly
    lxind = GetParameterIndex(model.name,'lx');
    eccind = GetParameterIndex(model.name,'ecc');   
    [ul_vec ul_ind] = unique([model.params(:,lxind) model.params(:,eccind)],'rows','stable');
    modelcode = 2; % planar
    for j = 1:length(ul_ind)        
        D = model.params(ul_ind(j),GetParameterIndex(model.name,'di')); % diffusivity
        a = model.params(ul_ind(j),GetParameterIndex(model.name,'lx'))/2;  % lx/2
        Ecc = model.params(ul_ind(j),GetParameterIndex(model.name,'ecc'));  % ecc        
        c = Ecc*a;
        roo=MMroots(dimA,a,modelcode); % x dir
        [protocol.S_x{j} protocol.A_x{j} ]=MatrixSA(qunit,a,dimA,roo,modelcode);       
        protocol.S_y{j} = protocol.S_x{j};
        protocol.A_y{j} = protocol.A_x{j};
        
        roo = MMroots(dimA,c,modelcode); % x dir
        [protocol.S_z{j} protocol.A_z{j} ]=MatrixSA(qunit,c,dimA,roo,modelcode);   
        % from roo MatrixR takes only the eigenvalues which are the same
        % no mater what length
        protocol.R_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
        protocol.R_y{j}=protocol.R_x{j};
        protocol.R_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,D);
        if protocol.diff==1
            % lx = ly pert
            apert=a*(1+protocol.pert);
            roo=MMroots(dimA,apert,modelcode);
            [protocol.Sperta_x{j} protocol.Aperta_x{j} ]=MatrixSA(qunit,apert,dimA,roo,modelcode);
            protocol.Rperta_x{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D);  
            protocol.Sperta_y{j} = protocol.Sperta_x{j};
            protocol.Aperta_y{j} = protocol.Aperta_x{j};
            protocol.Rperta_y{j} = protocol.Rperta_x{j};
             cpert = apert*Ecc;
             roo=MMroots(dimA,cpert,modelcode);
             [protocol.Sperta_z{j} protocol.Aperta_z{j} ]=MatrixSA(qunit,cpert,dimA,roo,modelcode);
             protocol.Rperta_z{j}=MatrixR(protocol.tau,cpert,dimA,roo,modelcode,D);  
          
            %Ecc pert
            Eccpert = Ecc*(1+protocol.pert);
            cpert=Eccpert*a;
            roo=MMroots(dimA,cpert,modelcode);
            [protocol.SpertEcc_z{j} protocol.ApertEcc_z{j} ]=MatrixSA(qunit,cpert,dimA,roo,modelcode);
            protocol.RpertEcc_z{j}=MatrixR(protocol.tau,cpert,dimA,roo,modelcode,D); 
         
            % Dpert
            Dpert=D*(1+protocol.pert);
            protocol.RpertD_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
            protocol.RpertD_y{j}=protocol.RpertD_x{j}; 
            protocol.RpertD_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,Dpert); 
        end
    end
end
 protocol.diff=0; %Used in Synth
 end



% numCombs = size(model.restriction_params,1);
% % Calculating the matrices
% numres = length(model.restriction_name); % model.restriction is a cell array of restricted compartments
% for j = 1:numCombs
%      params=model.restriction_params(j,:); % parameters associated to restriction;
%      %can be different to the total model parameters when there is a composite model
%     for nr = 1:numres
%     if strcmp(model.restriction_name{nr},'Cylinder') % two parameters - diffusivity and radius
%         modelcode = 1; % cylinder
%         D = params(1); % diffusivity
%         a = params(2); % radius
%         roo=MMroots(dimA,a,modelcode); % cylinder
%         [protocol.S_cyl{j} protocol.A_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,a,dimA,roo,modelcode);
%         protocol.R_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
%         if protocol.diff==1
%             % Radius apert
%             apert=a*(1+protocol.pert);
%             roo=MMroots(dimA,apert,modelcode);
%             [protocol.Sperta_cyl{j} protocol.Aperta_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,apert,dimA,roo,modelcode);
%             protocol.Rperta_cyl{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D); 
%             % Dpert
%             Dpert=D*(1+protocol.pert);
%             roo=MMroots(dimA,a,modelcode);
%             protocol.RpertD_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
%         end
%         
%     elseif strcmp(model.restriction_name{nr},'FiniteCylinder')
%         D = params(1); % diffusivity
%         a = params(2); % radius
%         l = params(3)/2; % half cylinder length; in the code the parallel planes are situated between -l and l
%         roo=MMroots(dimA,a,1); % cylinder
%         [protocol.S_cyl{j} protocol.A_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,a,dimA,roo,1);
%         protocol.R_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,1,D);
%         
%         roo=MMroots(dimA,l,2); % parallel planes
%         [protocol.S_plane{j} protocol.A_plane{j} ]=MatrixSA(qunit,l,dimA,roo,2);
%         protocol.R_plane{j}=MatrixR(protocol.tau,l,dimA,roo,2,D);
%          if protocol.diff==1
%             % Radius apert
%             apert=a*(1+protocol.pert);
%             roo=MMroots(dimA,apert,1);
%             [protocol.Sperta_cyl{j} protocol.Aperta_cyl{j} protocol.kDn_cyl{j}]=MatrixSA(qunit,apert,dimA,roo,1);
%             protocol.Rperta_cyl{j}=MatrixR(protocol.tau,apert,dimA,roo,1,D); 
%             % Dpert
%             % cylinder
%             Dpert=D*(1+protocol.pert);
%             protocol.RpertD_cyl{j}=MatrixR(protocol.tau,a,dimA,roo,1,Dpert); 
%             % plane
%             roo=MMroots(dimA,l,2);
%             protocol.RpertD_plane{j}=MatrixR(protocol.tau,l,dimA,roo,2,Dpert); 
%             % L perp
%             lpert=l*(1+protocol.pert);
%              roo=MMroots(dimA,lpert,2);
%             [protocol.Spertl_plane{j} protocol.Apertl_plane{j}]=MatrixSA(qunit,lpert,dimA,roo,2);
%             protocol.Rpertl_plane{j}=MatrixR(protocol.tau,lpert,dimA,roo,2,D); 
%          end        
%         
%     elseif strcmp(model.restriction_name{nr},'Sphere')
%         modelcode = 0; % sphere
%         D = params(1); % diffusivity
%         a = params(2); % radius
%         roo=MMroots(dimA,a,modelcode); % sphere
%         [Stemp Atemp protocol.kDn_sph]=MatrixSA(qunit,a,dimA,roo,modelcode);
%         protocol.A_sph{j,1} = Atemp{1}; % A0
%         protocol.A_sph{j,2} = Atemp{2}; % A90
%         protocol.S_sph{j,1} = Stemp{1}; % S0
%         protocol.S_sph{j,2} = Stemp{2}; % S90
%         protocol.R_sph{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
%         if protocol.diff==1
%             % Radius apert
%             apert=a*(1+protocol.pert);
%             roo=MMroots(dimA,apert,modelcode);
%             [Stemp Atemp ]=MatrixSA(qunit,apert,dimA,roo,modelcode);
%             protocol.Aperta_sph{j,1} = Atemp{1}; % A0
%             protocol.Aperta_sph{j,2} = Atemp{2}; % A90
%             protocol.Sperta_sph{j,1} = Stemp{1}; % S0
%             protocol.Sperta_sph{j,2} = Stemp{2}; % S90
%             protocol.Rperta_sph{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D); 
%             % Dpert
%             Dpert=D*(1+protocol.pert);
%             roo=MMroots(dimA,a,modelcode);
%             protocol.RpertD_sph{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
%         end
%     elseif strcmp(model.restriction_name{nr},'Cuboid')
%         modelcode = 2; % parallel planes E = Ex*Ey*Ez
%         D = params(1); % diffusivity
%         a = params(2)/2; % lx/2
%         b = params(3)/2; % ly/2
%         c = params(4)/2; % lz/2
%         roo=MMroots(dimA,a,modelcode); % x dir
%         [protocol.S_x{j} protocol.A_x{j} ]=MatrixSA(qunit,a,dimA,roo,modelcode);
%         roo = MMroots(dimA,b,modelcode); % x dir
%         [protocol.S_y{j} protocol.A_y{j} ]=MatrixSA(qunit,b,dimA,roo,modelcode);
%         roo = MMroots(dimA,c,modelcode); % x dir
%         [protocol.S_z{j} protocol.A_z{j} ]=MatrixSA(qunit,c,dimA,roo,modelcode);   
%         % from roo MatrixR takes only the eigenvalues which are the same
%         % no mater what length
%         protocol.R_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
%         protocol.R_y{j}=MatrixR(protocol.tau,b,dimA,roo,modelcode,D);
%         protocol.R_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,D);
%         if protocol.diff==1
%             % lx pert
%             apert=a*(1+protocol.pert);
%             roo=MMroots(dimA,apert,modelcode);
%             [protocol.Sperta_x{j} protocol.Aperta_x{j} ]=MatrixSA(qunit,apert,dimA,roo,modelcode);
%             % ly pert
%             bpert=b*(1+protocol.pert);
%             roo=MMroots(dimA,bpert,modelcode);
%             [protocol.Spertb_y{j} protocol.Apertb_y{j} ]=MatrixSA(qunit,bpert,dimA,roo,modelcode);
%             %lz pert
%             cpert=c*(1+protocol.pert);
%             roo=MMroots(dimA,cpert,modelcode);
%             [protocol.Spertc_z{j} protocol.Apertc_z{j} ]=MatrixSA(qunit,cpert,dimA,roo,modelcode);
%       
%             protocol.Rperta_x{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D); 
%             protocol.Rpertb_y{j}=MatrixR(protocol.tau,bpert,dimA,roo,modelcode,D); 
%             protocol.Rpertc_z{j}=MatrixR(protocol.tau,cpert,dimA,roo,modelcode,D); 
%             % Dpert
%             Dpert=D*(1+protocol.pert);
%             protocol.RpertD_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
%             protocol.RpertD_y{j}=MatrixR(protocol.tau,b,dimA,roo,modelcode,Dpert); 
%             protocol.RpertD_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,Dpert); 
%         end
%         
%     elseif strcmp(model.restriction_name{nr},'CuboidSym')  
%         modelcode = 2; % parallel planes E = Ex*Ey*Ez
%         D = params(1); % diffusivity
%         a = params(2)/2; % lx/2
%         Ecc = params(3);  % the last parameter is Ecc = lz/lx
%         c = Ecc*a;
%         roo=MMroots(dimA,a,modelcode); % x dir
%         [protocol.S_x{j} protocol.A_x{j} ]=MatrixSA(qunit,a,dimA,roo,modelcode);       
%         protocol.S_y{j} = protocol.S_x{j};
%         protocol.A_y{j} = protocol.A_x{j};
%         
%         roo = MMroots(dimA,c,modelcode); % x dir
%         [protocol.S_z{j} protocol.A_z{j} ]=MatrixSA(qunit,c,dimA,roo,modelcode);   
%         % from roo MatrixR takes only the eigenvalues which are the same
%         % no mater what length
%         protocol.R_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,D);
%         protocol.R_y{j}=protocol.R_x{j};
%         protocol.R_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,D);
%         if protocol.diff==1
%             % lx = ly pert
%             apert=a*(1+protocol.pert);
%             roo=MMroots(dimA,apert,modelcode);
%             [protocol.Sperta_x{j} protocol.Aperta_x{j} ]=MatrixSA(qunit,apert,dimA,roo,modelcode);
%             protocol.Rperta_x{j}=MatrixR(protocol.tau,apert,dimA,roo,modelcode,D);  
%             protocol.Sperta_y{j} = protocol.Sperta_x{j};
%             protocol.Aperta_y{j} = protocol.Aperta_x{j};
%             protocol.Rperta_y{j} = protocol.Rperta_x{j};
%              cpert = apert*Ecc;
%              roo=MMroots(dimA,cpert,modelcode);
%              [protocol.Sperta_z{j} protocol.Aperta_z{j} ]=MatrixSA(qunit,cpert,dimA,roo,modelcode);
%              protocol.Rperta_z{j}=MatrixR(protocol.tau,cpert,dimA,roo,modelcode,D);  
%           
%             %Ecc pert
%             Eccpert = Ecc*(1+protocol.pert);
%             cpert=Eccpert*a;
%             roo=MMroots(dimA,cpert,modelcode);
%             [protocol.SpertEcc_z{j} protocol.ApertEcc_z{j} ]=MatrixSA(qunit,cpert,dimA,roo,modelcode);
%             protocol.RpertEcc_z{j}=MatrixR(protocol.tau,cpert,dimA,roo,modelcode,D); 
%          
%             % Dpert
%             Dpert=D*(1+protocol.pert);
%             protocol.RpertD_x{j}=MatrixR(protocol.tau,a,dimA,roo,modelcode,Dpert); 
%             protocol.RpertD_y{j}=protocol.RpertD_x{j}; 
%             protocol.RpertD_z{j}=MatrixR(protocol.tau,c,dimA,roo,modelcode,Dpert); 
%         end
%     else
%         error('Unknown restriction for Matrix Method')
%     end
%     end
% end    
% protocol.diff=0; %Used in Synth
% end
