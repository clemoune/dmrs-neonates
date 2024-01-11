function [wave_form]= wave_form(protocol)
%
% camino.m--------------------------------------------------------------
% Discrete diffusion gradient waveform
% 
% [wave_form]= wave_form(protocol)
% 
% Description: Returns the discretized gradient waveforms of the diffusion 
% measurements specified in protoco that can be used as input to synthesize
% the diffusion signal using Matrix Method formalism
% Suported pulse sequences in protocol: PGSE, OGSE, SWOGSE, TWOGSE, SWOGSE_3D
% dPFG, STEAM, Helical
%
% Parameters:
% wave_form - the discretized gradient waveform, size [N, 3*K], 
%       where N is the number of measurements and K is the number of
%       gradient points in each direction.
%       [G11x G11y G11z G12x ... G1Kx G1Ky G1Kz]
%       ......
%       [GN1x GN1y GN1z GN2x ... GNKx GNKy GNKz]
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal. 
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

if ~isfield(protocol,'tau')
    protocol.tau = 1E-4;
end
if strcmp(protocol.pulseseq, 'PGSE') || strcmp(protocol.pulseseq, 'OGSE') || strcmp(protocol.pulseseq, 'SWOGSE') || strcmp(protocol.pulseseq, 'dSWOGSE') || ...
        strcmp(protocol.pulseseq, 'isoPGSE') ||strcmp(protocol.pulseseq, 'TWOGSE') || strcmp(protocol.pulseseq, 'SWOGSE_3D') % sequences based on spin echo sequences
    if isfield(protocol,'totalmeas')
        M = protocol.totalmeas;
    else
        M = length(protocol.smalldel);
    end
    total_time = max(ceil((protocol.smalldel+protocol.delta)/protocol.tau)).*protocol.tau+protocol.tau;
    time_vec = 0:protocol.tau:total_time;
    wav = zeros(3,length(time_vec));
    wave_form=zeros(M,length(wav(:)));
end

if(strcmp(protocol.pulseseq, 'PGSE'))
    for j = 1:M
        dG =protocol.G(j).*protocol.grad_dirs(j,:);
        for i=2:length(time_vec)
            if(time_vec(i)-protocol.smalldel(j)<-1E-10)
               wav(:,i) = dG;
               i_smalldel = i+1;
               i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
            elseif(time_vec(i)-protocol.delta(j)<-1E-10)
               wav(:,i) = 0;
               i_delta = i;
            elseif i-i_delta <=i_smalldel % operate on integers
                if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests                    
                    wav(:,i) = -wav(:,i-i_delta+1);
                else                    
                    wav(:,i) = -wav(:,i-i_delta);                    
                end
            else
                wav(:,i) = 0;
            end
        end
         wave_form(j,:) = wav(:);
    end
% elseif (strcmp(protocol.pulseseq, 'dPGSE'))
%     total_time = max(ceil((protocol.smalldel+2*protocol.delta+protocol.tm)/protocol.tau)).*protocol.tau+protocol.tau;
%     time_vec = 0:protocol.tau:total_time;
%  
%     M = size(protocol.grad_dirs1,1);
%     wave_form=zeros(M,length(time_vec)*3);
%     
%     
%     for j = 1:M
%         wav = zeros(3,length(time_vec));
%         dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:)';
%         dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:)';
%        
%         
%       index_smalldel = find(time_vec>0 & time_vec<protocol.smalldel(j));
%       wav(:,index_smalldel) = repmat(dG1,1,length(index_smalldel));
%       i_delta = find((time_vec - protocol.delta(j)) < protocol.tau+1E-10 & (time_vec - protocol.delta(j)) >0);
%       wav(:,i_delta+1:i_delta+length(index_smalldel)) = -repmat(dG1,1,length(index_smalldel));
%       
%       index_smalldel2= find((time_vec>protocol.delta(j)+ protocol.tm(j)) & (time_vec<protocol.smalldel(j)+protocol.delta(j)+protocol.tm(j)));
%       wav(:,index_smalldel2) = -repmat(dG2,1,length(index_smalldel2));
%       i_delta2 = find((time_vec - 2*protocol.delta(j)-protocol.tm(j)) < protocol.tau+1E-10 & (time_vec - 2*protocol.delta(j)-protocol.tm(j)) >0);
%       wav(:,i_delta2+1:i_delta2+length(index_smalldel2)) = repmat(dG2,1,length(index_smalldel2));  
%       
%   
%       wave_form(j,:) = wav(:);
%     end
    
elseif (strcmp(protocol.pulseseq, 'tPGSE'))
    total_time = max(ceil((protocol.smalldel+3*protocol.delta+2*protocol.tm)/protocol.tau)).*protocol.tau+3*protocol.tau;
    time_vec = 0:protocol.tau:total_time;
    
    M = size(protocol.grad_dirs1,1); 
    wave_form=zeros(M,length(time_vec)*3);
      
    
    for j = 1:M
        wav = zeros(3,length(time_vec));
        dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:)';
        dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:)';
        dG3 =protocol.G3(j).*protocol.grad_dirs3(j,:)';
        
      index_smalldel = find(time_vec>0 & time_vec<protocol.smalldel(j));
      wav(:,index_smalldel) = repmat(dG1,1,length(index_smalldel));
      i_delta = find((time_vec - protocol.delta(j)) < protocol.tau+1E-10 & (time_vec - protocol.delta(j)) >0);
      wav(:,i_delta+1:i_delta+length(index_smalldel)) = -repmat(dG1,1,length(index_smalldel));
      
      index_smalldel2= find((time_vec>protocol.delta(j)+ protocol.tm(j)) & (time_vec<protocol.smalldel(j)+protocol.delta(j)+protocol.tm(j)));
      wav(:,index_smalldel2) = -repmat(dG2,1,length(index_smalldel2));
      i_delta2 = find((time_vec - 2*protocol.delta(j)-protocol.tm(j)) < protocol.tau+1E-10 & (time_vec - 2*protocol.delta(j)-protocol.tm(j)) >0);
      wav(:,i_delta2+1:i_delta2+length(index_smalldel2)) = repmat(dG2,1,length(index_smalldel2));  
      
      index_smalldel3= find((time_vec>2*protocol.delta(j)+ 2*protocol.tm(j)) & (time_vec<protocol.smalldel(j)+2*protocol.delta(j)+2*protocol.tm(j)));
      wav(:,index_smalldel3) = repmat(dG3,1,length(index_smalldel3));
      i_delta3 = find((time_vec - 3*protocol.delta(j)-2*protocol.tm(j)) <protocol.tau+1E-10 & (time_vec - 3*protocol.delta(j)-2*protocol.tm(j)) >0);
      wav(:,i_delta3+1:i_delta3+length(index_smalldel3)) = -repmat(dG3,1,length(index_smalldel3));  
      
      wave_form(j,:) = wav(:);
    end
elseif(strcmp(protocol.pulseseq, 'isoPGSE'))
    for j = 1:M
        dG1 =protocol.G(j).*protocol.grad_dirs1(j,:);
        dG2 =protocol.G(j).*protocol.grad_dirs2(j,:);
        dG3 =protocol.G(j).*protocol.grad_dirs3(j,:);
        for i=2:length(time_vec)
            if(time_vec(i)-protocol.smalldel(j)/6<-1E-10)                
               wav(:,i) = dG1';
            elseif(time_vec(i)-2*protocol.smalldel(j)/6<-1E-10)   
               wav(:,i) = -dG1;
            elseif(time_vec(i)-3*protocol.smalldel(j)/6<-1E-10)   
               wav(:,i) = dG2;   
            elseif(time_vec(i)-4*protocol.smalldel(j)/6<-1E-10)   
               wav(:,i) = -dG2;
            elseif(time_vec(i)-5*protocol.smalldel(j)/6<-1E-10)   
               wav(:,i) = dG3;
            elseif(time_vec(i)-protocol.smalldel(j)<-1E-10)   
               wav(:,i) = -dG3;   
               i_smalldel = i+1;
               i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
            elseif(time_vec(i)-protocol.delta(j)<-1E-10)
               wav(:,i) = 0;
               i_delta = i;
            elseif i-i_delta <=i_smalldel % operate on integers
                if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests                    
                    wav(:,i) = -wav(:,i-i_delta+1);
                else                    
                    wav(:,i) = -wav(:,i-i_delta);                    
                end
            else
                wav(:,i) = 0;
            end
        end
         wave_form(j,:) = wav(:);
    end    
    

elseif(strcmp(protocol.pulseseq, 'SWOGSE')) 
     if(isfield(protocol,'omega'))
        for j = 1:M
        dG =protocol.G(j).*protocol.grad_dirs(j,:);
         
            if ~isfield(protocol,'phase')  % no phase            
                for i=2:length(time_vec)
                    if(time_vec(i)-protocol.smalldel(j)<-1E-10)       
                        wav(:,i) =dG*(-1).^floor(protocol.omega(j)./pi()*time_vec(i)-1E-10);
                        i_smalldel = i+1;
                        i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)-protocol.delta(j)<-1E-10)    
                        wav(:,i) = 0;
                        i_delta = i;
                   elseif i-i_delta <=i_smalldel % operate on integers
                       if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                             if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta+1);
                             else
                                 wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                             end
                        else
                             if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta);
                             else
                                wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                             end
                        end
                    else
                        wav(:,i) = 0;
                    end  
                end            
            else           
                for i=2:length(time_vec)                     
                    if(time_vec(i)-protocol.smalldel(j)<-1E-10)
                          wav(:,i) =  dG*(-1).^floor((protocol.omega(j)*time_vec(i)-protocol.phase(j))./pi()-1E-10);
                          i_smalldel = i+1;
                          i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)-protocol.delta(j)<-1E-10)
                          wav(:,i) = 0;
                          i_delta =  i;
                    elseif i-i_delta <=i_smalldel % operate on integers
                        if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                             if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta+1);
                             else
                                 wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                             end
                        else
                             if  ~isfield(protocol,'mirror') || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta);
                             else
                                wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                             end
                        end
                    else
                         wav(:,i) = 0;
                    end
                end           
            end
         wave_form(j,:) = wav(:);
        end
      else
        error('the protocol does not contain enough information to generate SWOGSE');
     end
elseif(strcmp(protocol.pulseseq, 'OGSE')) 
    if ~isfield(protocol,'apodisedcos') || protocol.apodisedcos == 0
     if(isfield(protocol,'omega'))
        for j = 1:M
        dG =protocol.G(j).*protocol.grad_dirs(j,:);
         
            if ~isfield(protocol,'phase')  % no phase            
                for i=2:length(time_vec)
                    if(time_vec(i)-protocol.smalldel(j)<-1E-10)       
                        wav(:,i) =dG*sin(protocol.omega(j).*time_vec(i));
                        i_smalldel = i+1;
                        i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)-protocol.delta(j)<-1E-10)    
                        wav(:,i) = 0;
                        i_delta = i;
                   elseif i-i_delta <=i_smalldel % operate on integers
                       if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                             if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta+1);
                             else
                                 wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                             end
                        else
                             if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta);
                             else
                                wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                             end
                        end
                    else
                        wav(:,i) = 0;
                    end  
                end            
            else           
                for i=2:length(time_vec)                     
                    if(time_vec(i)-protocol.smalldel(j)<-1E-10)
                          wav(:,i) =  dG*sin(protocol.omega(j)*time_vec(i)-protocol.phase(j));
                          i_smalldel = i+1;
                          i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)-protocol.delta(j)<-1E-10)
                          wav(:,i) = 0;
                          i_delta =  i;
                    elseif i-i_delta <=i_smalldel % operate on integers
                        if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                             if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta+1);
                             else
                                 wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                             end
                        else
                             if  ~isfield(protocol,'mirror') || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta);
                             else
                                wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                             end
                        end
                    else
                         wav(:,i) = 0;
                    end
                end           
            end
         wave_form(j,:) = wav(:);
        end
      else
        error('the protocol does not contain enough information to generate SWOGSE');
     end 
    else % returns the waveforms for apodised cosine (Does MRM 2003)
        
        for j = 1:M
        dG =protocol.G(j).*protocol.grad_dirs(j,:);
        N = floor(protocol.omega(j).*protocol.smalldel(j)/pi+0.0000000001); 
        if protocol.omega(j).*protocol.smalldel(j)/pi - floor(protocol.omega(j).*protocol.smalldel(j)/pi+0.0000000001) > 1E-4; 
            error('apodised cosine works only for waveforms with integer number of lobebs. To use this set the appropriate omega.')
        end
        if isfield(protocol,'phase')
            if abs(protocol.phase -pi/2) > 1E-4 || abs(protocol.phase + pi/2) > 1E-4 
                error('apodisation works only for cosine waveforms, i.e. phase = pi/2 or -pi/2; for a different phase protocol.apodisedcos must be 0')
            end
        end
        for i=1:length(time_vec)

            if(time_vec(i)-pi/protocol.omega(j)/2<-1E-10)
               wav(:,i) = dG*sin(2*protocol.omega(j)*time_vec(i));
            elseif(time_vec(i)-(protocol.smalldel(j)-pi/protocol.omega(j)/2)<-1E-10)
               wav(:,i) = -dG*sin(protocol.omega(j)*time_vec(i)-pi/2);
            elseif(time_vec(i)-protocol.smalldel(j)<-1E-10)
                wav(:,i) = (-1).^N.*dG*sin(2*protocol.omega(j)*(time_vec(i)-protocol.smalldel(j))+pi);
                i_smalldel = i+1;
                i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
            elseif(time_vec(i)-protocol.delta(j)<-1E-10)
                wav(:,i) = 0;
                 i_delta =  i;
           elseif i-i_delta <=i_smalldel % operate on integers
                if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                     if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                        wav(:,i) = -wav(:,i-i_delta+1);
                     else
                         wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                     end
                else
                     if  ~isfield(protocol,'mirror') || protocol.mirror ==0 % repeated
                        wav(:,i) = -wav(:,i-i_delta);
                     else
                        wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                     end
                end
            else
                 wav(:,i) = 0;
            end
        end
         wave_form(j,:) = wav(:);
        end        
    end
     
elseif strcmp(protocol.pulseseq,'TWOGSE')  % does not have phase yet
    if ~isfield(protocol,'apodisedcos') || protocol.apodisedcos == 0
    if(isfield(protocol,'slew_rate') && isfield(protocol,'omega') )
        rt_vec = protocol.G./protocol.slew_rate; % rise time  
        for j = 1:M  
            if ~isfield(protocol,'phase')  % no phase  
                rt = rt_vec(j);
                Nt = floor(protocol.smalldel(j).*protocol.omega(j)./pi+0.00000000001);
                dG =protocol.G(j).*protocol.grad_dirs(j,:);
                
                    for i=2:length(time_vec)    
                        if(time_vec(i)<Nt.*pi/protocol.omega(j))
                            it = floor(protocol.omega(j)*time_vec(i)./pi()+0.00000001);
                            if( time_vec(i)<it *pi/protocol.omega(j)+rt)
                                wav(:,i) = dG./rt.*(-1).^it.*(time_vec(i)-it.*pi./protocol.omega(j));

                            elseif( time_vec(i)<(it+1)*pi/protocol.omega(j)-rt)
                                wav(:,i) = dG*(-1)^it;
                            else                   
                                wav(:,i) = dG./rt.*(-1).^it.*((it+1).*pi./protocol.omega(j)-time_vec(i));
                            end
                             i_smalldel = i;
                             i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                        elseif(time_vec(i)<protocol.smalldel(j)) 
                            if (protocol.smalldel(j)-Nt.*pi/protocol.omega(j)>2*rt)  % last oscillation longer than 2 *rt
                                if( time_vec(i)<Nt.*pi/protocol.omega(j)+rt)
                                    wav(:,i) = dG./rt.*(-1).^Nt.*(time_vec(i)-Nt.*pi./protocol.omega(j));

                                elseif( time_vec(i)<protocol.smalldel(j)-rt)
                                    wav(:,i) = dG*(-1)^Nt;
                                else
                                    wav(:,i) = dG./rt.*(-1).^Nt.*(protocol.smalldel(j)-time_vec(i));
                                end
                            else  % last oscillation shorter than 2 *rt
                                if( time_vec(i)<Nt.*pi/protocol.omega(j)+(protocol.smalldel(j)-Nt.*pi/protocol.omega(j))/2)
                                    wav(:,i) = dG./rt.*(-1).^Nt.*(time_vec(i)-Nt.*pi./protocol.omega(j));  
                                else
                                    wav(:,i) = dG./rt.*(-1).^Nt.*(protocol.smalldel(j)-time_vec(i));
                                end
                            end
                            i_smalldel = i;
                            i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                        elseif(time_vec(i)<protocol.delta(j))
                            wav(:,i) = 0;
                            i_delta = i;
                        elseif i-i_delta <=i_smalldel % operate on integers
                            if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                                 if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                    wav(:,i) = -wav(:,i-i_delta+1);
                                 else
                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                                 end
                            else
                                 if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                    wav(:,i) = -wav(:,i-i_delta);
                                 else
                                    wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                                 end
                             end                        

                        else
                            wav(:,i) = 0;
                        end
                    end
              
            else % waveforms with phase
                rt = rt_vec(j);                
                dG =protocol.G(j).*protocol.grad_dirs(j,:);
                if(protocol.phase(j)>= 0)
                    phdelay = protocol.phase(j)./protocol.omega(j); % time delay due to phase shift
                else
                    phdelay = ((pi()-abs(protocol.phase(j)))./protocol.omega(j));
                    dG = -dG;
                end
                Nt = floor((protocol.smalldel(j)-phdelay).*protocol.omega(j)./pi+0.00000000001);
                    for i=2:length(time_vec)  
                        if time_vec(i)<phdelay
                            if phdelay > 2* rt % duration of the phase shift lrger than 2*rt
                                if(time_vec(i)-rt<-1E-10)
                                     wav(:,i) = dG./rt.*time_vec(i);
                                elseif(time_vec(i)-(phdelay-rt)<-1E-10)
                                      wav(:,i) = dG;
                                else
                                      wav(:,i) = dG./rt.*(phdelay-time_vec(i));
                                end
                            else
                                if(time_vec(i)< phdelay/2)
                                    wav(:,i) = dG./rt.*time_vec(i);                               
                                else
                                    wav(:,i) = dG./rt.*(phdelay-time_vec(i));
                                end
                            end
                            i_smalldel = i;
                            i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                            
                        elseif(time_vec(i)<Nt.*pi/protocol.omega(j)+ phdelay)
                            it = floor(protocol.omega(j)*(time_vec(i)-phdelay)./pi()+0.00000001);
                            if( time_vec(i)<it *pi/protocol.omega(j)+rt+phdelay)
                                wav(:,i) = dG./rt.*(-1).^(it+1).*(time_vec(i)-it.*pi./protocol.omega(j)-phdelay);

                            elseif( time_vec(i)<(it+1)*pi/protocol.omega(j)-rt +phdelay)
                                wav(:,i) = dG*(-1)^(it+1);
                            else                   
                                wav(:,i) = dG./rt.*(-1).^(it+1).*((it+1).*pi./protocol.omega(j)+phdelay-time_vec(i));
                            end
                             i_smalldel = i;
                            i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                        elseif(time_vec(i)<protocol.smalldel(j)) 
                            if (protocol.smalldel(j)-Nt.*pi/protocol.omega(j)-phdelay>2*rt) % last oscillation longer than 2 *rt
                                if( time_vec(i)<Nt.*pi/protocol.omega(j)+rt+ phdelay)
                                    wav(:,i) = dG./rt.*(-1).^(Nt+1).*(time_vec(i)-Nt.*pi./protocol.omega(j)-phdelay);

                                elseif( time_vec(i)<protocol.smalldel(j)-rt)
                                    wav(:,i) = dG*(-1)^(Nt+1);
                                else
                                    wav(:,i) = dG./rt.*(-1).^(Nt+1).*(protocol.smalldel(j)-time_vec(i));
                                end
                            else  % last oscillation shorter than 2 *rt
                                if( time_vec(i)<protocol.smalldel(j) -(protocol.smalldel(j)-Nt.*pi/protocol.omega(j)-phdelay)/2)
                                    wav(:,i) = dG./rt.*(-1).^(Nt+1).*(time_vec(i)-Nt.*pi./protocol.omega(j)-phdelay);  
                                else
                                    wav(:,i) = dG./rt.*(-1).^(Nt+1).*(protocol.smalldel(j)-time_vec(i));
                                end
                            end
                            i_smalldel = i;
                            i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                        elseif(time_vec(i)<protocol.delta(j))
                            wav(:,i) = 0;
                            i_delta = i;
                        elseif i-i_delta <=i_smalldel % operate on integers
                            if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                                 if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                    wav(:,i) = -wav(:,i-i_delta+1);
                                 else
                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                                 end
                            else
                                 if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                    wav(:,i) = -wav(:,i-i_delta);
                                 else
                                    wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                                 end
                             end                        

                        else
                            wav(:,i) = 0;
                        end
                    end
                
            end
            wave_form(j,:) = wav(:);
            
        end
    else
        error('the protocol does not contain enough information to generate TWOGSE');
    end 
    else % apodised trapezoid; waveforms described in Ianus et al JMR 2012. The area of the first lobe is half the area of a full lobe;
        % to achieve this, the duration of the waveform is increased to
        % smalldel + rt
        rt_vec = protocol.G./protocol.slew_rate; % rise time  
        for j = 1:M
            rt = rt_vec(j);
            dG =protocol.G(j).*protocol.grad_dirs(j,:);
            if protocol.omega(j).*protocol.smalldel(j)/pi - floor(protocol.omega(j).*protocol.smalldel(j)/pi+0.0000000001) > 1E-4; 
                error('apodised cosine works only for waveforms with integer number of lobebs. To use this option set the appropriate omega.')
            end
            if isfield(protocol,'phase')
                if abs(protocol.phase -pi/2) > 1E-4 & abs(protocol.phase + pi/2) > 1E-4 
                    error('apodisation works only for cosine waveforms, i.e. phase = pi/2 or -pi/2; for a different phase protocol.apodisedcos must be 0')
                end
            end
             for i=1:length(time_vec)
                if(time_vec(i)-(pi/protocol.omega(j)/2+rt/2)<-1E-10)
                    if(time_vec(i)-rt<-1E-10)
                         wav(:,i) = dG./rt.*time_vec(i);
                    elseif(time_vec(i)-(pi/2/protocol.omega(j)+rt/2-rt)<-1E-10)
                          wav(:,i) = dG;
                    else
                          wav(:,i) = dG./rt.*(pi/2./protocol.omega(j)+rt/2-time_vec(i));
                    end
                elseif(time_vec(i)-(protocol.smalldel(j)-pi/protocol.omega(j)/2+rt/2)<-1E-10)
                    it = floor(protocol.omega(j)*(time_vec(i)-pi/protocol.omega(j)/2-rt/2)./pi+0.00000000001);
                    if( (time_vec(i)-pi/protocol.omega(j)/2-rt/2)-(it *pi/protocol.omega(j)+rt)<-1E-10)
                        wav(:,i) = -dG./rt.*(-1).^it.*((time_vec(i)-pi/protocol.omega(j)/2-rt/2)-it.*pi./protocol.omega(j));
                    elseif((time_vec(i)-pi/protocol.omega(j)/2-rt/2)-((it+1)*pi/protocol.omega(j)-rt)<-1E-10)
                        wav(:,i) = -dG*(-1)^it; 
                    else
                        wav(:,i) = -dG./rt.*(-1).^it.*((it+1).*pi./protocol.omega(j)-(time_vec(i)-pi/protocol.omega(j)/2-rt/2));
                    end
                elseif(time_vec(i)-(protocol.smalldel(j)+rt)<-1E-10)
                     it = floor(protocol.omega(j)*(time_vec(i)-pi/protocol.omega(j)/2-rt/2)./pi+0.00000000001);
                     if( (time_vec(i)-pi/protocol.omega(j)/2-rt/2)-(it *pi/protocol.omega(j)+rt)<-1E-10)
                        wav(:,i) = -dG./rt.*(-1).^it.*((time_vec(i)-pi/protocol.omega(j)/2-rt/2)-it.*pi./protocol.omega(j));
                    elseif(time_vec(i) -protocol.smalldel(j)<-1E-10)
                        wav(:,i) = -dG*(-1)^it; 
                    else
                        wav(:,i) = -dG./rt.*(-1).^it.*(protocol.smalldel(j)+rt-time_vec(i));
                     end
                    i_smalldel = i;
                    i_delta = i;
                elseif(time_vec(i)-protocol.delta(j)+rt<-1E-10)
                    wav(:,i) = 0;
                    i_delta = i;
                elseif i-i_delta <=i_smalldel % operate on integers
                    if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta+1);
                         else
                             wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                         end
                    else
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end
                     end                        

                else
                    wav(:,i) = 0;
                end
             end
             wave_form(j,:) = wav(:);
        end
        
    end
elseif strcmp(protocol.pulseseq,'dSWOGSE')    
    if isfield(protocol,'grad_dirs1') && isfield(protocol,'grad_dirs2') 
        total_time = max(ceil((protocol.smalldel+2*protocol.delta+protocol.tm)/protocol.tau)).*protocol.tau+protocol.tau;
        time_vec = 0:protocol.tau:total_time;
        wav = zeros(3,length(time_vec));
        wave_form = zeros(M,length(wav(:))); 
        for j = 1:M
        dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
        dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);
        i_delta2 = 0;
         i_smalldel2 = 0;
         i_delta1 = 0;
         i_smalldel1 = 0; 
            if ~isfield(protocol,'phase')  % no phase            
                for i=2:length(time_vec)
                    if(time_vec(i)<protocol.smalldel(j))       
                        wav(:,i) =dG1*(-1).^floor((protocol.omega(j)*time_vec(i))./pi()-1E-10);
                        i_smalldel1 = i+1;
                        i_delta1 = i_smalldel1; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)<protocol.delta(j))  
                         wav(:,i) = 0;           
                        i_delta1 = i; % for the case when delta-smalldel = 0; needed  for some tests        

                    elseif i-i_delta1 <=i_smalldel1 % operate on integers
                       if i_delta1 == i_smalldel1 % needed when delta-smalldel = 0; needed  for some tests                           
                              if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta1+1);
                             else
                                wav(:,i) = -wav(:,i_smalldel1-(i-i_delta1));
                             end                          
                       else                             
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta1);
                         else
                            wav(:,i) = -wav(:,i_smalldel1-(i-i_delta1)+1);
                         end                         
                       end
                       i_tm = i;
                    elseif time_vec(i)<protocol.tm(j)+protocol.delta(j)+1E-10
                         wav(:,i) = 0;   
                         i_tm = i;
                    elseif(time_vec(i)<protocol.smalldel(j)+protocol.tm(j)+protocol.delta(j))       
                        wav(:,i) =-dG2*(-1).^floor((protocol.omega(j)*(time_vec(i)-protocol.tm(j)-protocol.delta(j)))./pi-1E-10);
                        i_smalldel2 = i+1;
                        i_delta2 = i_smalldel2; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)<2*protocol.delta(j)+protocol.tm(j))  
                         wav(:,i) = 0;           
                        i_delta2 = i; % for the case when delta-smalldel = 0; needed  for some tests        
            %       
                    elseif i-i_delta2 <=i_smalldel2-i_tm % operate on integers
                       if i_delta2 == i_smalldel2 % needed when delta-smalldel = 0; needed  for some tests      
                            if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta2+1+i_tm);
                             else
                                wav(:,i) = -wav(:,i_smalldel2-(i-i_delta2+i_tm));
                             end 

                       else 
                           if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                             wav(:,i) = -wav(:,i-i_delta2+i_tm);
                           else
                                wav(:,i) = -wav(:,i_smalldel2-(i-i_delta2+i_tm)+1);
                           end                

                       end 

                    else
                       wav(:,i) = 0;
                    end 
                end            
            else           
                for i=2:length(time_vec)                     
                    if(time_vec(i)<protocol.smalldel(j))       
                        wav(:,i) =dG1*(-1).^floor((protocol.omega(j)*time_vec(i)-protocol.phase(j))./pi()-1E-10);
                        i_smalldel1 = i+1;
                        i_delta1 = i_smalldel1; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)<protocol.delta(j))  
                         wav(:,i) = 0;           
                        i_delta1 = i; % for the case when delta-smalldel = 0; needed  for some tests        

                    elseif i-i_delta1 <=i_smalldel1 % operate on integers
                       if i_delta1 == i_smalldel1 % needed when delta-smalldel = 0; needed  for some tests                           
                              if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta1+1);
                             else
                                wav(:,i) = -wav(:,i_smalldel1-(i-i_delta1));
                             end                          
                       else                             
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta1);
                         else
                            wav(:,i) = -wav(:,i_smalldel1-(i-i_delta1)+1);
                         end                         
                       end
                       i_tm = i;
                    elseif time_vec(i)<protocol.tm(j)+protocol.delta(j)+1E-10
                         wav(:,i) = 0;   
                         i_tm = i;
                    elseif(time_vec(i)<protocol.smalldel(j)+protocol.tm(j)+protocol.delta(j))       
                        wav(:,i) =-dG2*(-1).^floor((protocol.omega(j)*(time_vec(i)-protocol.tm(j)-protocol.delta(j))-protocol.phase(j))./pi-1E-10);
                        i_smalldel2 = i+1;
                        i_delta2 = i_smalldel2; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(time_vec(i)<2*protocol.delta(j)+protocol.tm(j))  
                         wav(:,i) = 0;           
                        i_delta2 = i; % for the case when delta-smalldel = 0; needed  for some tests        
            %       
                    elseif i-i_delta2 <=i_smalldel2-i_tm % operate on integers
                       if i_delta2 == i_smalldel2 % needed when delta-smalldel = 0; needed  for some tests      
                            if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                                wav(:,i) = -wav(:,i-i_delta2+1+i_tm);
                             else
                                wav(:,i) = -wav(:,i_smalldel2-(i-i_delta2+i_tm));
                             end 

                       else 
                           if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                             wav(:,i) = -wav(:,i-i_delta2+i_tm);
                           else
                                wav(:,i) = -wav(:,i_smalldel2-(i-i_delta2+i_tm)+1);
                           end                
 M = size(protocol.grad_dirs1,1);
                       end 

                    else
                       wav(:,i) = 0;
                    end 
                end           
            end
         wave_form(j,:) = wav(:);
        end
      else
        error('the protocol does not contain enough information to generate dSWOGSE');
    end
elseif strcmp(protocol.pulseseq,'DODE')    
    if isfield(protocol,'grad_dirs1') && isfield(protocol,'grad_dirs2') 
       if ~isfield(protocol,'slew_rate')
        if ~isfield(protocol,'phase') || protocol.phase == 0           
            
            M = size(protocol.grad_dirs1,1);
            total_time = max(ceil((protocol.smalldel+protocol.tm+protocol.smalldel)/protocol.tau)).*protocol.tau+protocol.tau;
            time_vec = 0:protocol.tau:total_time;
           
            wave_form = zeros(M,length(wav(:))); 
            np1 = floor(protocol.smalldel./2./protocol.Nosc./protocol.tau);
            np2 = floor(protocol.smalldel./2./protocol.Nosc./protocol.tau);
            for j = 1:M
                 wav = zeros(3,length(time_vec));
            dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
            dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

            wf1 = repmat([1 -1],np1(j),1);
            wf1 = repmat(wf1(:)',1,protocol.Nosc(j));
            wf2 = repmat([-1 1],np2(j),1);
            wf2 = repmat(wf2(:)',1,protocol.Nosc(j));
            int = zeros(3,floor(protocol.tm(j)./protocol.tau));

            wav(:,2:1+(length(wf1)+length(int)+length(wf2))) = [repmat(dG1',1,length(wf1)).*repmat(wf1,3,1) ...
                int repmat(dG2',1,length(wf2)).*repmat(wf2,3,1)];

             wave_form(j,:) = wav(:);
            end
        elseif protocol.phase == pi/2 || protocol.phase == -pi/2
             M = size(protocol.grad_dirs1,1);
            total_time = max(ceil((protocol.smalldel+protocol.tm+protocol.smalldel)/protocol.tau)).*protocol.tau+protocol.tau;
            time_vec = 0:protocol.tau:total_time;
        
            wave_form = zeros(M,3*length(time_vec)); 
            np1 = floor(protocol.smalldel./4./protocol.Nosc./protocol.tau);
            np2 = floor(protocol.smalldel./4./protocol.Nosc./protocol.tau);
            for j = 1:M
                 wav = zeros(3,length(time_vec));
            dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
            dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

            wf1 = repmat([1 -1 -1 1],np1(j),1);
            wf1 = repmat(wf1(:)',1,floor(protocol.Nosc(j)));
            
            if rem(protocol.Nosc(j),1) 
                wf1end = repmat([1 -1],np1(j),1);
                wf1 = [wf1 wf1end(:)'];
            end
            wf2 = repmat([-1 1 1 -1],np2(j),1);
            wf2 = repmat(wf2(:)',1,floor(protocol.Nosc(j)));
            if rem(protocol.Nosc(j),1) 
                wf2end = repmat([-1 1],np2(j),1);
                wf2 = [wf2 wf2end(:)'];
            end
            int = zeros(3,floor(protocol.tm(j)./protocol.tau));

            wav(:,2:1+(length(wf1)+length(int)+length(wf2))) = [repmat(dG1',1,length(wf1)).*repmat(wf1,3,1) ...
                int repmat(dG2',1,length(wf2)).*repmat(wf2,3,1)];

             wave_form(j,:) = wav(:);
            end
            
        else
            error('for DODE sequences the phase must be either 0 or pi/2')
        end
       else
           if ~isfield(protocol,'phase') || protocol.phase == 0           
            rt1 = protocol.G1./protocol.slew_rate;
            rt2 = protocol.G2./protocol.slew_rate;
            M = size(protocol.grad_dirs1,1);
            total_time = max(ceil((protocol.smalldel+protocol.tm+protocol.smalldel)./protocol.tau)).*protocol.tau+protocol.tau;
            time_vec = 0:protocol.tau:total_time;
           
            wave_form = zeros(M,3*length(time_vec)); 
            np1 = floor(protocol.smalldel./2./protocol.Nosc./protocol.tau);
            np2 = floor(protocol.smalldel./2./protocol.Nosc./protocol.tau);
            np1_rt = floor(rt1./protocol.tau);
            np2_rt = floor(rt2./protocol.tau);
            np1 = np1 - 2*np1_rt;
            np2 = np2 - 2*np2_rt;            
                    

            
            for j = 1:M
                 wav = zeros(3,length(time_vec));
           wf1up = linspace(0,1,np1_rt(j)+1);
           wf1up = wf1up(2:end);
           wf1down = linspace(1,0,np1_rt(j)+1);
           wf1down = wf1down(2:end);

           wf2up = linspace(0,1,np2_rt(j)+1);
           wf2up = wf2up(2:end);
           wf2down = linspace(1,0,np2_rt(j)+1);
           wf2down = wf2down(2:end);    
                
            dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
            dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

            wf1 = [wf1up ones(1,np1(j)) wf1down -wf1up -ones(1,np1(j)) -wf1down];
            wf1 = repmat(wf1,1,protocol.Nosc(j));
            wf2 = [-wf2up -ones(1,np2(j)) -wf2down wf2up ones(1,np2(j)) wf2down];
            wf2 = repmat(wf2,1,protocol.Nosc(j));
            int = zeros(3,floor(protocol.tm(j)./protocol.tau));

            wav(:,2:1+(length(wf1)+length(int)+length(wf2))) = [repmat(dG1',1,length(wf1)).*repmat(wf1,3,1) ...
                int repmat(dG2',1,length(wf2)).*repmat(wf2,3,1)];

             wave_form(j,:) = wav(:);
            end
           elseif protocol.phase == pi/2 || protocol.phase == -pi/2
            rt1 = protocol.G1./protocol.slew_rate;
            rt2 = protocol.G2./protocol.slew_rate;
            M = size(protocol.grad_dirs1,1);
            total_time = max(ceil((protocol.smalldel+protocol.tm+protocol.smalldel+ rt1 + rt2)./protocol.tau)).*protocol.tau+protocol.tau;
            time_vec = 0:protocol.tau:total_time;
            
            wave_form = zeros(M,3*length(time_vec)); 
            
            
            np11 = floor(((protocol.smalldel./4./protocol.Nosc)+rt1/2)./protocol.tau);
            np12 = floor((protocol.smalldel./2./protocol.Nosc)./protocol.tau);
            np21 = floor(((protocol.smalldel./4./protocol.Nosc)+rt2/2)./protocol.tau);
            np22 = floor((protocol.smalldel./2./protocol.Nosc)./protocol.tau);
            
            np1_rt = floor(rt1./protocol.tau);
            np2_rt = floor(rt2./protocol.tau);
            np11 = np11 - 2*np1_rt;
            np12 = np12 - 2*np1_rt;
            np21 = np21 - 2*np2_rt;
            np22 = np22 - 2*np2_rt;  
                    
            Nhp = (protocol.Nosc*2); % number of half periods
            
            for j = 1:M
           wav = zeros(3,length(time_vec));
           
           wf1up = linspace(0,1,np1_rt(j)+1);
           wf1up = wf1up(2:end);
           wf1down = linspace(1,0,np1_rt(j)+1);
           wf1down = wf1down(2:end);

           wf2up = linspace(0,1,np2_rt(j)+1);
           wf2up = wf2up(2:end);
           wf2down = linspace(1,0,np2_rt(j)+1);
           wf2down = wf2down(2:end);    
                
            dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
            dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

              if mod(Nhp(j),2) == 0 % even number of half periods            
            
            wf10 = [-wf1up -ones(1,np12(j)) -wf1down wf1up ones(1,np12(j)) wf1down];
            wf1 = [wf1up ones(1,np11(j)) wf1down repmat(wf10,1,(protocol.Nosc(j)-1)) ...
                -wf1up -ones(1,np12(j)) -wf1down wf1up ones(1,np11(j)) wf1down ];
            wf20 = [wf2up ones(1,np22(j)) wf2down -wf2up -ones(1,np22(j)) -wf1down];
            wf2 = [-wf2up -ones(1,np21(j)) -wf2down repmat(wf20,1,(protocol.Nosc(j)-1)) wf2up ones(1,np22(j)) wf2down -wf2up -ones(1,np21(j)) -wf2down ];
            int = zeros(3,floor(protocol.tm(j)./protocol.tau));

                       
              else
            wf10 = [-wf1up -ones(1,np12(j)) -wf1down wf1up ones(1,np12(j)) wf1down];
            wf1 = [wf1up ones(1,np11(j)) wf1down repmat(wf10,1,(floor(protocol.Nosc(j)))) -wf1up -ones(1,np11(j)) -wf1down ];
            wf20 = [wf2up ones(1,np22(j)) wf2down -wf2up -ones(1,np22(j)) -wf1down];
            wf2 = [-wf2up -ones(1,np21(j)) -wf2down repmat(wf20,1,(floor(protocol.Nosc(j)))) wf2up ones(1,np21(j)) wf2down ];
            int = zeros(3,floor(protocol.tm(j)./protocol.tau));

              end
               wav(:,2:1+(length(wf1)+length(int)+length(wf2))) = [repmat(dG1',1,length(wf1)).*repmat(wf1,3,1) ...
                int repmat(dG2',1,length(wf2)).*repmat(wf2,3,1)];
              wave_form(j,:) = wav(:);
           end
           
           end
       end
      else
        error('the protocol does not contain enough information to generate DODE');
     end     
elseif strcmp(protocol.pulseseq,'DDE') || strcmp(protocol.pulseseq,'dPGSE') ||  strcmp(protocol.pulseseq,'dPFG')       
    if isfield(protocol,'grad_dirs1') && isfield(protocol,'grad_dirs2') 
        
        if ~isfield(protocol,'slew_rate')
        
            if sum(protocol.tm >= protocol.smalldel) == length(protocol.tm)
            M = size(protocol.grad_dirs1,1);
            total_time = max(ceil((2*protocol.delta+protocol.tm+protocol.smalldel)/protocol.tau)).*protocol.tau+10*protocol.tau;
            time_vec = 0:protocol.tau:total_time;
           
            wave_form = zeros(M,3*length(time_vec)); 
            np1 = floor(protocol.smalldel./protocol.tau);
            np2 = floor(protocol.smalldel./protocol.tau);
            for j = 1:M
              wav = zeros(3,length(time_vec));   
            dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
            dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

            wf1 = ones(3,np1(j));
            wf2 = ones(3,np2(j));
            int1 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);
            int2 = zeros(3,floor((protocol.tm(j)-protocol.smalldel(j))./protocol.tau)+1);
            int3 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);

            wav(:,2:1+(2*size(wf1,2)+size(int1,2)+2*size(wf2,2)+size(int2,2)+size(int3,2))) = ...
                [repmat(dG1',1,size(wf1,2)).*wf1 int1 -repmat(dG1',1,size(wf1,2)).*wf1 ...
                int2 -repmat(dG2',1,size(wf2,2)).*wf2 int3 repmat(dG2',1,size(wf2,2)).*wf2];

             wave_form(j,:) = wav(:);
            end
            elseif sum(protocol.tm == 0) == length(protocol.tm)

                 M = size(protocol.grad_dirs1,1);
            total_time = max(ceil((2*protocol.delta+protocol.tm+protocol.smalldel)/protocol.tau)).*protocol.tau+10*protocol.tau;
            time_vec = 0:protocol.tau:total_time;
            
            wave_form = zeros(M,3*length(time_vec)); 
            np1 = floor(protocol.smalldel./protocol.tau);
            np2 = floor(protocol.smalldel./protocol.tau);
            for j = 1:M
            wav = zeros(3,length(time_vec));    
            dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
            dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

            wf1 = ones(3,np1(j));
            wf2 = ones(3,np2(j));
            int1 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);

            int3 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);

            wav(:,2:1+(2*size(wf1,2)+size(int1,2)+size(wf2,2)+size(int3,2))) = ...
                [repmat(dG1',1,size(wf1,2)).*wf1 int1 -repmat(dG1'+dG2',1,size(wf1,2)).*wf1 ...
                  int3 repmat(dG2',1,size(wf2,2)).*wf2];

             wave_form(j,:) = wav(:);
            end

            else error ('waveform not yet defined for 0< tm <smalldel')
            end
        else
            if sum(protocol.tm >= protocol.smalldel) == length(protocol.tm)
                rt1 = protocol.G1./protocol.slew_rate;
                rt2 = protocol.G2./protocol.slew_rate;           
                
                M = size(protocol.grad_dirs1,1);
                total_time = max(ceil((2*protocol.delta+protocol.tm+protocol.smalldel)/protocol.tau)).*protocol.tau+10*protocol.tau;
                time_vec = 0:protocol.tau:total_time;
               
                wave_form = zeros(M,3*length(time_vec)); 
                np1 = floor(protocol.smalldel./protocol.tau);
                np2 = floor(protocol.smalldel./protocol.tau);
                np1_rt = floor(rt1./protocol.tau);
                np2_rt = floor(rt2./protocol.tau);
                np1 = np1 - 2*np1_rt;
                np2 = np2 - 2*np2_rt;   
                
                
                for j = 1:M
                     wav = zeros(3,length(time_vec));
                   wf1up = linspace(0,1,np1_rt(j)+1);
                   wf1up = wf1up(2:end);
                   wf1down = linspace(1,0,np1_rt(j)+1);
                   wf1down = wf1down(2:end);

                   wf2up = linspace(0,1,np2_rt(j)+1);
                   wf2up = wf2up(2:end);
                   wf2down = linspace(1,0,np2_rt(j)+1);
                   wf2down = wf2down(2:end);    
                    
                    dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
                    dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);

                    wf1 = [wf1up ones(1,np1(j)) wf1down];
                    wf2 = [wf2up ones(1,np2(j)) wf2down];
                    int1 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);
                    int2 = zeros(3,floor((protocol.tm(j)-protocol.smalldel(j))./protocol.tau)+1);
                    int3 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);

                    wav(:,2:1+(2*length(wf1)+size(int1,2)+2*length(wf2)+size(int2,2)+size(int3,2))) = ...
                        [repmat(dG1',1,length(wf1)).*repmat(wf1,3,1) int1 -repmat(dG1',1,length(wf1)).*repmat(wf1,3,1) ...
                        int2 -repmat(dG2',1,length(wf2)).*repmat(wf2,3,1) int3 repmat(dG2',1,length(wf2)).*repmat(wf2,3,1)];

                     wave_form(j,:) = wav(:);
                end
            elseif sum(protocol.tm == 0) == length(protocol.tm)
                 rt1 = protocol.G1./protocol.slew_rate;
                rt2 = protocol.G2./protocol.slew_rate;     
                 M = size(protocol.grad_dirs1,1);
                total_time = max(ceil((2*protocol.delta+protocol.tm+protocol.smalldel)/protocol.tau)).*protocol.tau+10*protocol.tau;
                time_vec = 0:protocol.tau:total_time;
              
                wave_form = zeros(M,3*length(time_vec)); 
               np1 = floor(protocol.smalldel./protocol.tau);
                np2 = floor(protocol.smalldel./protocol.tau);
                np1_rt = floor(rt1./protocol.tau);
                np2_rt = floor(rt2./protocol.tau);
                np1 = np1 - 2*np1_rt;
                np2 = np2 - 2*np2_rt;                  
                
                
                for j = 1:M
                     wav = zeros(3,length(time_vec)); 
                    dG1 =protocol.G1(j).*protocol.grad_dirs1(j,:);
                    dG2 =protocol.G2(j).*protocol.grad_dirs2(j,:);                   
                   
                    rt12 = norm(dG1'+dG2')./protocol.slew_rate;
                    np12_rt = floor(rt12./protocol.tau);
                    np12 = np1 - 2*np12_rt;
                    
                    wf1up = linspace(0,1,np1_rt(j)+1);
                   wf1up = wf1up(2:end);
                   wf1down = linspace(1,0,np1_rt(j)+1);
                   wf1down = wf1down(2:end);

                   wf2up = linspace(0,1,np2_rt(j)+1);
                   wf2up = wf2up(2:end);
                   wf2down = linspace(1,0,np2_rt(j)+1);
                   wf2down = wf2down(2:end);    
                   
                   wf12up = linspace(0,1,np12_rt+1);
                   wf12up = wf12up(2:end);
                   wf12down = linspace(1,0,np12_rt+1);
                   wf12down = wf12down(2:end);


                    wf1 = [wf1up ones(1,np1(j)) wf1down];
                    wf2 = [wf2up ones(1,np2(j)) wf2down];
                    wf12 = [wf12up ones(1,np12) wf12down];
      
                    int1 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);

                    int3 = zeros(3,floor((protocol.delta(j)-protocol.smalldel(j))./protocol.tau)+1);

                    wav(:,2:1+(2*length(wf1)+length(int1)+length(wf2)+length(int3))) = ...
                        [repmat(dG1',1,length(wf1)).*wf1 int1 -repmat(dG1'+dG2',1,length(wf12)).*wf12 ...
                          int3 repmat(dG2',1,length(wf2)).*wf2];

                     wave_form(j,:) = wav(:);
                end

            else error ('waveform not yet defined for 0< tm <smalldel')
            end
            
        end
        
        
      else
        error('the protocol does not contain enough information to generate DDE');
     end       
elseif strcmp(protocol.pulseseq,'SWOGSE_3D')
    tau=protocol.tau;
    K = floor((protocol.smalldel+1E-10)./tau)+1;
    dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
    wave_form=zeros(M,3*(max(K)*2+max(dstmp)));
    if ~isfield(protocol,'angle') || protocol.angle == 4
         Gx = protocol.Gx';
         Gy = protocol.Gy';
         Gz = protocol.Gz';
         if ~isfield(protocol,'phix') && ~isfield(protocol,'phiy') && ~isfield(protocol,'phiz')
            phix = zeros(1,M); phiy = zeros(1,M); phiz = zeros(1,M);
        else
            phix = protocol.phix; phiy = protocol.phiy;  phiz = protocol.phiz;
         end
        for m=1:M              

                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-phiy(m))./pi-1E-10);
                sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-phiz(m))./pi-1E-10);
                vec1 = zeros(size(sign_vec1)); vec2 = zeros(size(sign_vec2));  vec3 = zeros(size(sign_vec3)); 
                vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end)=0;
                vec2(sign_vec2>0) = Gy(m); vec2(sign_vec2<0) = -Gy(m); vec2(1)=0; vec2(end)=0;
                vec3(sign_vec3>0) = Gz(m); vec3(sign_vec3<0) = -Gz(m); vec3(1)=0; vec3(end)=0;

                  if ~isfield(protocol,'mirror') || protocol.mirror == 0
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)]; 
                    Gy_vec = [vec2 zeros(1,dstmp(m)) -(vec2)]; 
                    Gz_vec = [vec3 zeros(1,dstmp(m)) -(vec3)]; 
                  else
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)]; 
                    Gy_vec = [vec2 zeros(1,dstmp(m)) -fliplr(vec2)]; 
                    Gz_vec = [vec3 zeros(1,dstmp(m)) -fliplr(vec3)]; 
                   
                  end                       
                wave_form(m,1:3:length(Gx_vec)*3) = Gx_vec;
                wave_form(m,2:3:length(Gy_vec)*3) = Gy_vec;
                wave_form(m,3:3:length(Gz_vec)*3) = Gz_vec;

        end
    elseif protocol.angle == 1
         Gx = protocol.Gx';
         if ~isfield(protocol,'phix')
            phix = zeros(1,M);
        else
            phix = protocol.phix;
         end
         for m = 1:M
            time_vec = tau*(0:K(m)-1);
            sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
            vec1 = zeros(size(sign_vec1));
            vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end)=0;
               if ~isfield(protocol,'mirror') || protocol.mirror == 0
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)];                    
               else
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)]; 
               end                       
             wave_form(m,1:3:length(Gx_vec)*3) = Gx_vec;
         end

    elseif protocol.angle == 2
         Gx = protocol.Gx';
      Gy = protocol.Gy';
       
        if ~isfield(protocol,'phix') && ~isfield(protocol,'phiy')
            phix = zeros(1,M); phiy = zeros(1,M);
        else
            phix = protocol.phix; phiy = protocol.phiy; 
        end

            for m=1:M                 

                    time_vec = tau*(0:K(m)-1);
                    sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                    sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-phiy(m))./pi-1E-10);
                    vec1 = zeros(size(sign_vec1)); vec2 = zeros(size(sign_vec2)); 
                    vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end)=0;
                    vec2(sign_vec2>0) = Gy(m); vec2(sign_vec2<0) = -Gy(m); vec2(1)=0; vec2(end)=0;


                   if ~isfield(protocol,'mirror') || protocol.mirror == 0
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)]; 
                    Gy_vec = [vec2 zeros(1,dstmp(m)) -(vec2)]; 
                   
                  else
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)]; 
                    Gy_vec = [vec2 zeros(1,dstmp(m)) -fliplr(vec2)]; 
                   
                  end                       
                wave_form(m,1:3:length(Gx_vec)*3) = Gx_vec;
                wave_form(m,2:3:length(Gy_vec)*3) = Gy_vec;
            end
    elseif protocol.angle == 3
         Gx = protocol.Gx';
          Gz = protocol.Gz';
       
         if ~isfield(protocol,'phix') && ~isfield(protocol,'phiz')
            phix = zeros(1,M); phiz = zeros(1,M);
        else
            phix = protocol.phix; phiz = protocol.phiz; 
        end

            for m=1:M                 

                    time_vec = tau*(0:K(m)-1);
                    sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                    sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-phiz(m))./pi-1E-10);
                    vec1 = zeros(size(sign_vec1)); vec2 = zeros(size(sign_vec2)); 
                    vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end) = 0;
                    vec2(sign_vec2>0) = Gz(m); vec2(sign_vec2<0) = -Gz(m); vec2(1)=0; vec2(end) = 0;


                   if ~isfield(protocol,'mirror') || protocol.mirror == 0
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)];                    
                    Gz_vec = [vec2 zeros(1,dstmp(m)) -(vec2)]; 
                  else
                    Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)];                   
                    Gz_vec = [vec2 zeros(1,dstmp(m)) -fliplr(vec2)]; 
                   
                  end                       
                wave_form(m,1:3:length(Gx_vec)*3) = Gx_vec;
                wave_form(m,3:3:length(Gz_vec)*3) = Gz_vec;
                   
            end
    else error('Unknown protocol.angle')
    end
        
%     total_time = max(ceil((protocol.smalldel+protocol.delta)/protocol.tau)).*protocol.tau+protocol.tau;
%     time_vec = 0:protocol.tau:total_time;
%     wav = zeros(3,length(time_vec));
%     wave_form=zeros(M,length(wav(:)));
%     if ~isfield(protocol,'angle') || protocol.angle == 4
%         if(isfield(protocol,'omegax') && isfield(protocol,'omegay') && isfield(protocol,'omegaz'))
%             if (~isfield(protocol,'phix') && ~isfield(protocol,'phiz') && ~isfield(protocol,'phiy'))
%                 for j = 1:M
%                 dG =[protocol.Gx(j); protocol.Gy(j); protocol.Gz(j)];
%                    
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([protocol.omegax(j)*time_vec(i)./pi; protocol.omegay(j)*time_vec(i)./pi();...
%                                 protocol.omegaz(j)*time_vec(i)./pi()]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                  if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end   
%                    
%                    wave_form(j,:) = wav(:);   
%                 end
%             else
%                 for j = 1:M
%                 dG =[protocol.Gx(j); protocol.Gy(j); protocol.Gz(j)];
%                    
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([(protocol.omegax(j)*time_vec(i)-protocol.phix(j))./pi; (protocol.omegay(j)*time_vec(i)-protocol.phiy(j))./pi();...
%                                 (protocol.omegaz(j)*time_vec(i)-protocol.phiz(j))./pi()]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end 
%                     wave_form(j,:) = wav(:);
%                   
%                 end
%             end
%         else error('the protocol does not contain enough information to generate the waveform')
%         end
%     elseif protocol.angle == 1
%         if(isfield(protocol,'omegax'))
%             if (~isfield(protocol,'phix'))
%                 for j = 1:M
%                 dG =[protocol.Gx(j); 0; 0];
%                   
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([protocol.omegax(j)*time_vec(i)./pi; 0; 0]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                              if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end 
%                     wave_form(j,:) = wav(:);
%                 end
%                   
%             else
%                 for j = 1:M
%                 dG =[protocol.Gx(j); 0; 0];
%                   
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([(protocol.omegax(j)*time_vec(i)-protocol.phix(j))./pi; 0; 0]-1E-10);
%                             i_smalldel = i+1;
%                             i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                   else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                                 
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end 
%                    wave_form(j,:) = wav(:); 
%                   
%                 end
%             end
%         else error('the protocol does not contain enough information to generate the waveform');
%         end
%     elseif protocol.angle == 2
%         if(isfield(protocol,'omegax') && isfield(protocol,'omegay') )
%             if (~isfield(protocol,'phix')  && ~isfield(protocol,'phiy'))
%                 for j = 1:M
%                 dG =[protocol.Gx(j); protocol.Gy(j); 0];                   
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([protocol.omegax(j)*time_vec(i)./pi; protocol.omegay(j)*time_vec(i)./pi; 0]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end   
%                     wave_form(j,:) = wav(:);
%                 end               
%             else
%                 for j = 1:M
%                 dG =[protocol.Gx(j); protocol.Gy(j); 0];
%                    
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([(protocol.omegax(j)*time_vec(i)-protocol.phix(j))./pi; (protocol.omegay(j)*time_vec(i)-protocol.phiy(j))./pi; 0]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end   
%                     wave_form(j,:) = wav(:);
%                 end
%             end
%         else error('the protocol does not contain enough information to generate the waveform');
%         end
%     elseif protocol.angle == 3
%         if(isfield(protocol,'omegax')  && isfield(protocol,'omegaz'))
%             if (~isfield(protocol,'phix') && ~isfield(protocol,'phiz') )
%                 for j = 1:M
%                 dG =[protocol.Gx(j); 0; protocol.Gz(j)];
%                    
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([protocol.omegax(j)*time_vec(i)./pi; 0;...
%                                 protocol.omegaz(j)*time_vec(i)./pi]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end   
%                    
%                    wave_form(j,:) = wav(:);   
%                 end
%             else
%                 for j = 1:M
%                 dG =[protocol.Gx(j); 0; protocol.Gz(j)];
%                    
%                    for i=2:length(time_vec)
%                         if(time_vec(i)-protocol.smalldel(j)<-1E-10)  
%                             wav(:,i) = dG.*(-1).^floor([(protocol.omegax(j)*time_vec(i)-protocol.phix(j))./pi; 0;...
%                                 (protocol.omegaz(j)*time_vec(i)-protocol.phiz(j))./pi]-1E-10);
%                             i_smalldel = i+1;
%                         elseif (time_vec(i)-protocol.delta(j)<-1E-10) 
%                             wav(:,i) = 0;
%                             i_delta = i;
%                         elseif i-i_delta <=i_smalldel % operate on integers
%                             if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta+1);
%                                  else
%                                      wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
%                                  end
%                             else
%                                   if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
%                                     wav(:,i) = -wav(:,i-i_delta);
%                                  else
%                                     wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
%                                  end
%                             end
%                         else
%                             wav(:,i) = 0;
%                         end                     
%                    end   
%                    wave_form(j,:) = wav(:);
%                 end
%             end
%         else error('the protocol does not contain enough information to generate the waveform')
%         end
%     else error('Unknown protocol.angle')
%     end
elseif strcmp(protocol.pulseseq,'GEN') % creating a waveform when GENGx is specified   
  
    wave_form = protocol.G; 
   
elseif strcmp(protocol.pulseseq,'STEAM')    
    % stimulated echo sequence with imaging gradients along z direction and
    % diffusion gradients along grad_dirs
    smalldel = protocol.smalldel; % diffusion gradient duration
    sdelc = protocol.sdelc; % duration of crusher gradient
    sdels = protocol.sdels; % duration of slice select gradient
    tau = protocol.tau; % time interval for waveform discretisation
    G = protocol.G; % amplitude of the diffusion gradient
    Gc = protocol.Gc; % amplitude of the crusher gradient
    Gs = protocol.Gs; % amplitude of the slice gradient
    tau1 = protocol.tau1; % time interval between diffusion gradient and first crusher gradient
    tau2 = protocol.tau2; % time interval between the second crusher gradient and the 2nd diffusion gradient
    tm = protocol.tm; % mixing time between the 2 90deg rf pulses
    
    total_time = max(ceil((2*(smalldel+2*sdels+sdelc)+tau1+tau2+tm)/tau)).*tau+10*tau;
    M = length(protocol.smalldel);
    time_vec = 0:tau:total_time;   
    wave_form=zeros(M,3*length(time_vec));       

     for j = 1:M
        Nsdelc = round(sdelc(j)./tau);
        Nsdels = round(sdels(j)./tau);
        Ntau1 = round(tau1(j)./tau);
        Ntaum = round(tm(j)./tau);
        Ntau2 = round(tau2(j)./tau);
         wav = zeros(3,length(time_vec));
        dG =G(j).*protocol.grad_dirs(j,:);
        Nsmalldel = round(smalldel(j)./tau); % number of point in smalldel
        wav(:,2:1+Nsmalldel) = repmat(dG',1,Nsmalldel);
        wav(3,2+Nsmalldel+Ntau1:1+Nsmalldel+Ntau1+Nsdelc) = Gc(j);
        wav(3,2+Nsmalldel+Ntau1+Nsdelc:1+Nsmalldel+Ntau1+Nsdelc+Nsdels) = Gs(j);
        wav(3,2+Nsmalldel+Ntau1+Nsdelc+Nsdels+Ntaum:1+Nsmalldel+Ntau1+Nsdelc+2*Nsdels+Ntaum) = Gs(j);
        wav(3,2+Nsmalldel+Ntau1+Nsdelc+2*Nsdels+Ntaum:1+Nsmalldel+Ntau1+2*Nsdelc+2*Nsdels+Ntaum) = -Gc(j);
        wav(:,2+Nsmalldel+Ntau1+2*Nsdelc+2*Nsdels+Ntaum+Ntau2:1+2*Nsmalldel+Ntau1+2*Nsdelc+2*Nsdels+Ntaum+Ntau2) = -repmat(dG',1,Nsmalldel);       
       
        wave_form(j,:) = wav(:)';
     end    
    
% elseif strcmp(protocol.pulseseq,'dPFG') || strcmp(protocol.pulseseq,'dPGSE')
%     M = length(protocol.G);
%     tau = protocol.tau;
%     total_time = max(ceil((protocol.smalldel+protocol.delta+protocol.tm)/protocol.tau)).*protocol.tau+10*protocol.tau;
%     time_vec = 0:tau:total_time;   
%     wave_form=zeros(M,3*length(time_vec));
%     for j = 1:M
%     G = protocol.G(j);
%     smalldel = protocol.smalldel(j);
%     delta = protocol.delta(j);
%     tmix = protocol.tm(j);
%     theta = protocol.theta(j); 
%     phi = protocol.phi(j); % angles in spherical coord describing the orientation of the first pair of gradients in the lab coordinate system
%     theta1 = protocol.theta1(j);
%     phi1 = protocol.phi1(j);% angle in spherical coord describing the orientation of the second pair of gradients in the lab coordinate system
%    
% 
%     t1=smalldel; %time during first lobe of G1
%     t2=delta-smalldel; %between the end of the first lobe and the begining of the second lobe of G1
%     if tmix<smalldel
%       t3=tmix; % during just the second lobe of G1
%       t4=smalldel-tmix; % overlap of second lobe of G1 and first lobe of G2
%       t5=0; % between the end of the second lobe of G1 and begining of the first lobe of G2
%       t6=tmix;% just during the first lobe of G2
%     else
%       t3=smalldel;
%       t4=0;
%       t5=tmix-smalldel;
%       t6=smalldel;
%     end
%     t7=delta-smalldel;
%     t8=smalldel;
%     
%     tt1=round(t1/tau);
%     tt2=round(t2/tau);
%     tt3=round(t3/tau);
%     tt4=round(t4/tau);
%     tt5=round(t5/tau);
%     tt6=round(t6/tau);
%     tt7=round(t7/tau);
%     tt8=round(t8/tau);     
% 
%     G1v=G*[ [0 0 0]' ...
%         repmat([cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)]',1,tt1)  ...
%         repmat([0 0 0]',1,tt2)  ...
%         repmat(-[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)]',1,tt3)  ...
%         repmat(-[cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)]',1,tt4)  ...
%         repmat([0 0 0]',1,tt5)  ...
%         repmat([0 0 0]',1,tt6)  ...
%         repmat([0 0 0]',1,tt7)  ...
%         repmat([0 0 0]',1,tt8)  ...
%         [0 0 0]'];
%     
%     G2v=G*[ [ 0 0 0]' ...
%         repmat([0 0 0]',1,tt1)  ...
%         repmat([0 0 0]',1,tt2)  ...
%         repmat([0 0 0]',1,tt3)  ...
%         repmat([-cos(phi1)*sin(theta1) -sin(phi1)*sin(theta1) -cos(theta1)]',1,tt4)  ...
%         repmat([0 0 0]',1,tt5)  ...
%         repmat([-cos(phi1)*sin(theta1) -sin(phi1)*sin(theta1) -cos(theta1)]',1,tt6)  ...
%         repmat([0 0 0]',1,tt7)  ...
%         repmat([cos(phi1)*sin(theta1) sin(phi1)*sin(theta1) cos(theta1)]',1,tt8)  ...
%         [0 0 0]'];
%     
%     wav=G1v+G2v; 
%      wave_form(j,1:3*size(wav,2))=wav(:)';
%     end
elseif strcmp(protocol.pulseseq,'Helical')
    M = length(protocol.G);
    tau = protocol.tau;
    total_time = max(ceil((protocol.smalldel+protocol.delta)/protocol.tau)).*protocol.tau+10*protocol.tau;
    time_vec = 0:tau:total_time;   
    wave_form=zeros(M,3*length(time_vec));
    for j = 1:M
        G = protocol.G(j);       
        for i=2:length(time_vec)
            if(time_vec(i)-protocol.smalldel(j)<-1E-10)       
                wav(1,i) =G*cos(protocol.omega(j).*time_vec(i));
                wav(2,i) =G*sin(protocol.omega(j).*time_vec(i));
                wav(3,i) = G*time_vec(i).*protocol.slopez(j);
                i_smalldel = i+1;
                i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
            elseif(time_vec(i)-protocol.delta(j)<-1E-10)    
                wav(:,i) = 0;
                i_delta = i;
           elseif i-i_delta <=i_smalldel % operate on integers
               if i_delta == i_smalldel % needed when delta-smalldel = 0; needed  for some tests
                     if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                        wav(:,i) = -wav(:,i-i_delta+1);
                     else
                         wav(:,i) = -wav(:,i_smalldel-(i-i_delta));
                     end
                else
                     if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                        wav(:,i) = -wav(:,i-i_delta);
                     else
                        wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                     end
                end
            else
                wav(:,i) = 0;
            end  
        end
        wave_form(j,:) = wav(:)';
    end
   
    
else
    error('Unknown protocol.pulseseq');
end
